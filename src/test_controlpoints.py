import re
import sys

import matplotlib.pyplot as plt
import numpy as np


def read_points_from_off(filename: str):
    xs = []
    ys = []
    zs = []

    r = re.compile(r'([0-9.e-]+) ([0-9.e-]+) ([0-9.e-]+)')

    with open(filename, 'r') as file:
        for line in file.readlines():
            line = line.strip()

            if not line:
                if xs:
                    yield xs, ys, zs
                    xs = []
                    ys = []
                    zs = []
                continue

            m = r.match(line)
            if m:
                xs.append(float(m.group(1)))
                ys.append(float(m.group(2)))
                zs.append(float(m.group(3)))
            else:
                print(f'Cannot match line "{line}". Ignoring it.')


def next_color(color=None):
    if color == 'r':
        return 'g'
    if color == 'g':
        return 'b'
    if color == 'b':
        return 'y'
    else:
        return 'r'


def main(filename: str):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    all_xs = []
    all_ys = []
    all_zs = []

    c = None
    for xs, ys, zs in read_points_from_off(filename):
        all_xs.extend(xs)
        all_ys.extend(ys)
        all_zs.extend(zs)
        c = next_color(c)
        ax.scatter(xs, ys, zs, c=c, marker='o')
    ax.set_box_aspect((np.ptp(all_xs), np.ptp(all_ys), np.ptp(all_zs)))

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main('surfaces.off')