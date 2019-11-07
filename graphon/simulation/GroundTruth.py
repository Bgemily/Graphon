import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from scipy.stats import gamma, expon, norm


p = 30
Max_time = 100
alpha = .25
lamda = .1


# Generate node positions, specify cell types
np.random.seed(831)
data = np.random.rand(30, 2)
# cell_types = np.repeat([[0], [1]], [3, p-3], axis=0)
cell_types = np.repeat([[0], [1], [2]], [3, 4, p-7], axis=0)
data = np.concatenate((data, cell_types), axis=1)



# Generate edge times
T_mat = np.ones((p, p)) * -1


seed = 0
for i in range(p):
    for j in range(i+1, p):
        if cell_types[i]==cell_types[j]:
            edge_time = expon.rvs(scale=1/lamda, random_state=seed)
        elif 0 in {int(cell_types[i]), int(cell_types[j])}:
            edge_time = gamma.rvs(a=alpha, scale=1/lamda, random_state=seed)
        else: #cell types = {1,2}
            edge_time = norm.rvs(1, 1, random_state=seed)
        # edge_time = gamma.rvs(a=alpha, scale=1/lamda, random_state=seed) if cell_types[i] != cell_types[j] else expon.rvs(scale=1/lamda, random_state=seed)
        if edge_time < Max_time:
            T_mat[i, j] = edge_time
        T_mat[j, i] = T_mat[i, j]
        seed = seed + 1



if __name__ == '__main__': # visualization

    # nodes = (data[0,:2], data[1,:2], data[2, :2], data[3:,:2])
    # colors = ("red", "blue", "green", "gray")
    # groups = ("MN1", "MN2", "MN3", "Others")
    # markers = ("X", "X", "X", ".")
    nodes = (data[:3,:2], data[3:7,:2], data[7:,:2])
    colors = ("red", "blue", "gray")
    groups = ("Type1", "Type2", "Others")
    markers = ("X", "X", ".")


    for t in [.01, .1, 1, 10, 100]:

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        for node, color, group, marker in zip(nodes, colors, groups, markers):
            try:
                x, y = node[:,0], node[:,1]
            except:
                x, y = node
            ax.scatter(x, y, alpha=0.8, c=color, marker = marker, edgecolors='none', s=70, label=group)


        Is, Js = np.where(T_mat < t)

        for i, j in zip(Is, Js):
            line = lines.Line2D(data[[i,j],0], data[[i,j],1])
            line.set_linewidth(.4)
            # if 0 in {i,j}:
            #     line.set_color(colors[0])
            # elif 1 in {i,j}:
            #     line.set_color(colors[1])
            # elif 2 in {i,j}:
            #     line.set_color(colors[2])
            # else:
            #     line.set_color(colors[3])
            if 0 in {int(cell_types[i]),int(cell_types[j])}:
                line.set_color(colors[0])
            elif 1 in {int(cell_types[i]),int(cell_types[j])}:
                line.set_color(colors[1])
            else:
                line.set_color(colors[2])


            ax.add_line(line, )

        plt.title('Network(t={})'.format(t))
        plt.legend()
        plt.show()