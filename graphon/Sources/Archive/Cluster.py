import numpy as np
from scipy.spatial.distance import pdist
from scipy.stats import norm
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt




def intensity_est(poin_proc, h): # kernel: normal pdf
    n = len(poin_proc)
    return lambda t: 1/(n*h)*sum([norm.pdf((t-ti)/h) for ti in poin_proc])


def bandwidth_sele(poin_proc):
    time_range = max(poin_proc)-min(poin_proc)
    log_likes = []
    for h in np.arange(time_range/(2*len(poin_proc)),(time_range)/2, time_range/(2*len(poin_proc))):
        f_hat = intensity_est(poin_proc, h)
        log_like = sum(poin_proc*np.log(f_hat(poin_proc))+(1-poin_proc)*np.log(1-f_hat(poin_proc)))
        log_likes.append([log_like, h])

    log_likes = np.asarray(log_likes)

    h_opt = log_likes[np.argmax(log_likes[:,0]),1]
    return h_opt



def cro_corr(f1, f2):
    x = np.arange(-10, 100, 0.01)
    try:
        return np.linalg.norm(f1(x)-f2(x))*0.1
    except TypeError:
        f1, f2 = f1[0], f2[0]
        return np.linalg.norm(f1(x) - f2(x)) * 0.1


def clustering(T_mat):
    f_hats = []
    for row in T_mat:
        h_opt = bandwidth_sele(row)
        f_hats.append(intensity_est(row, h_opt))


    pdists = np.ones_like(T_mat)
    for i in range(len(f_hats)):
        for j in range(i+1, len(f_hats)):
            pdists[i,j] = cro_corr(f_hats[i], f_hats[j])
            pdists[j,i] = pdists[i,j]

    clusters = AgglomerativeClustering(n_clusters=3, affinity="precomputed", linkage="average").fit(pdists)

    return f_hats, clusters








if __name__ == '__main__':

    from GroundTruth import T_mat, data
    f_hats, clusters = clustering(T_mat)

    print(clusters.labels_)


    x = np.arange(0,100,0.01)

    # color_map1 = {'0': 'red', '1': 'black'}
    # color_map2 = {'0': 'pink', '1': 'lightgray'}

    color_map1 = {'0': 'red', '1': 'blue', '2': 'black'}
    color_map2 = {'0': 'pink', '1': 'lightblue', '2': 'lightgray'}

    # for type in [0,1]:
    #     nodes = np.where(data[:,2]==type)[0]
    #     mean_curve = np.mean([f_hats[node](x) for node in nodes], axis=0)
    #     sd_curve = np.std([f_hats[node](x) for node in nodes], axis=0)
    #     cl1, cl2 = color_map1[str(type)], color_map2[str(type)]
    #     label = "MNs" if type==0 else "Others"
    #     plt.plot(x, mean_curve, color=cl1, linestyle='dashed', linewidth=1, label=label)
    #     plt.fill_between(x, mean_curve-sd_curve, mean_curve+sd_curve, color=cl2)
    for type in [0,1,2]:
        nodes = np.where(data[:,2]==type)[0]
        mean_curve = np.mean([f_hats[node](x) for node in nodes], axis=0)
        sd_curve = np.std([f_hats[node](x) for node in nodes], axis=0)
        cl1, cl2 = color_map1[str(type)], color_map2[str(type)]
        label = "Type{}".format(type) if type != 2 else "Others"
        plt.plot(x, mean_curve, color=cl1, linestyle='dashed', linewidth=1, label=label)
        plt.fill_between(x, mean_curve-sd_curve, mean_curve+sd_curve, color=cl2)

    plt.title('Fitted intensity function')
    plt.legend()
    plt.xlabel('Time')
    plt.show()

