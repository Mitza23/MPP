import numpy as np
import matplotlib.pyplot as plt

file_path = './project/output/parallel_data_step_100.txt'

def read_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    dimensions = lines[0].strip().split()
    nx, ny, nz = int(dimensions[0]), int(dimensions[1]), int(dimensions[2])

    u = np.zeros((nx, ny, nz))
    v = np.zeros((nx, ny, nz))
    t = np.zeros((nx, ny, nz))

    for line in lines[1:]:
        data = line.strip().split()
        x, y, z = int(data[0]), int(data[1]), int(data[2])
        u_val, v_val, t_val = float(data[3]), float(data[4]), float(data[5])
        u[x, y, z] = u_val
        v[x, y, z] = v_val
        t[x, y, z] = t_val

    return nx, ny, nz, u, v, t

def plot_3d_surface(nx, ny, nz, data, title='Temperature'):
    X, Y, Z = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz))
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    kw = {
        'vmin': data.min(),
        'vmax': data.max(),
        'levels': np.linspace(data.min(), data.max(), 10),
    }
    # Plot contour surfaces
    _ = ax.contourf(
        X[:, :, 0], Y[:, :, 0], data[:, :, 0],
        zdir='z', offset=nz-1, **kw
    )
    _ = ax.contourf(
        X[0, :, :], data[0, :, :], Z[0, :, :],
        zdir='y', offset=0, **kw
    )
    C = ax.contourf(
        data[:, 0, :], Y[:, 0, :], Z[:, 0, :],
        zdir='x', offset=nx-1, **kw
    )

    xmin, xmax = X.min(), X.max()
    ymin, ymax = Y.min(), Y.max()
    zmin, zmax = Z.min(), Z.max()
    ax.set(xlim=[xmin, xmax], ylim=[ymin, ymax], zlim=[zmin, zmax])

    edges_kw = dict(color='0.4', linewidth=1, zorder=1e3)
    ax.plot([xmax, xmax], [ymin, ymax], 0, **edges_kw)
    ax.plot([xmin, xmax], [ymin, ymin], 0, **edges_kw)
    ax.plot([xmax, xmax], [ymin, ymin], [zmin, zmax], **edges_kw)

    ax.set(
        xlabel='X',
        ylabel='Y',
        zlabel='Z',
        zticks=[zmin, zmax // 2, zmax],
    )

    ax.view_init(40, -30)
    ax.set_box_aspect(None, zoom=0.9)

    fig.colorbar(C, ax=ax, fraction=0.02, pad=0.1, label=title)

    plt.show()

if __name__ == '__main__':
    nx, ny, nz, u, v, t = read_data(file_path)
    plot_3d_surface(nx, ny, nz, t, title='Temperature')
