import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # registers something for 3d plots


def plot_curve(curve, grid=True, separate=False, subdivision=0, output_file=None):
    """
    Draw a plot for peano multifractal 2d or 3d curve

    Args:
        grid:  draw grid?
        separate:  separate proto+pnums and subdivision protolines
        subdivision:  number of subdivision
        output_file:  subj
    """

    # render config
    pattern_names = 'abcd'
    colors = ['tab:cyan','orange','green','pink']
    box_aspect = (1, 1, 0.9)
    sub_linewidth = 0.8
    dpi = 300
    connect_linestyle = (0, (1, 1))  # dotted
    if separate:
        proto_style = {'linewidth': 2, 'alpha': 1}
    else:
        proto_style = {'linewidth': 5, 'alpha': 0.4}

    # patterns are stacked along x-axis
    dim = curve.dim
    assert dim <= 3
    div = curve.div
    kws = {} if dim == 2 else {'projection': '3d'}

    #fig = plt.figure(figsize=plt.figaspect(0.5))
    fig = plt.figure()
    axs = []
    ticks = [j/div for j in range(1, div)]
    print(ticks)
    for pnum in range(curve.pcount):
        pnum_axs = []
        if separate:
            ax1 = fig.add_subplot(curve.pcount, 2, pnum + 1, **kws)
            ax2 = fig.add_subplot(curve.pcount, 2, pnum + 1 + curve.pcount, **kws)
            pnum_axs = [ax1, ax2]
        else:
            pnum_axs = [fig.add_subplot(curve.pcount, 1, pnum + 1, **kws)]
        for idx, ax in enumerate(pnum_axs):
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            if grid:
                ax.set_xticks(ticks)
                ax.set_yticks(ticks)
            else:
                ax.set_xticks([])
                ax.set_yticks([])
            if curve.pcount > 1:
                if (not separate) or (separate and idx == 0):
                    ax.set_title('pattern ({})'.format(pattern_names[pnum]))
            if dim == 2:
                ax.set_aspect(1)
            elif dim == 3:
                ax.set_box_aspect(box_aspect)
                ax.set_zlim(0, 1)
                ax.set_zticklabels([])
                if grid:
                    ax.set_zticks(ticks)
                else:
                    ax.set_zticks([])
            if grid:
                ax.grid(linestyle='dotted')
            axs.append(ax)

    fig.subplots_adjust(wspace=-0.0,hspace=0.0)

    if subdivision:
        subcurve = curve.get_subdivision(subdivision)
    else:
        subcurve = curve

    for pnum, pattern in enumerate(curve.patterns):
        if not separate:
            glob_ax = detail_ax = axs[pnum]
        else:
            glob_ax = axs[2*pnum]
            detail_ax = axs[2*pnum+1]

        prev_line = None
        prev_color = None
        for cnum, (cube, spec) in enumerate(zip(pattern.proto, pattern.specs)):
            pname = pattern_names[spec.pnum]
            if spec.base_map.time_rev:
                pname = r'$\overline{\mathrm{' + pname + '}}$'
            if curve.pcount > 1:
                glob_ax.text(*[(cj + 0.37)/curve.div for cj in cube], s=pname, fontsize='medium', alpha=0.7)
            line = _get_proto_center_line(spec.base_map * subcurve.patterns[spec.pnum].proto)
            abs_line = [_get_abs_point(point, shift=cube, scale=1/curve.div) for point in line]
            if prev_line:
                med = tuple((aj+bj)/2 for aj, bj in zip(prev_line[-1], abs_line[0]))
                draw_lines(detail_ax, [prev_line[-1], med], linewidth=sub_linewidth, linestyle=connect_linestyle, color=prev_color)
                draw_lines(detail_ax, [med, abs_line[0]], linewidth=sub_linewidth, linestyle=connect_linestyle, color=colors[spec.pnum])
            draw_lines(detail_ax, abs_line, arrow=(cnum == curve.genus-1), linewidth=sub_linewidth, color=colors[spec.pnum])
            prev_line = abs_line
            prev_color = colors[spec.pnum]

        proto_col = colors[pnum] if separate else 'gray'
        draw_lines(glob_ax, _get_proto_center_line(pattern.proto), arrow=True, color=proto_col, **proto_style)

    if output_file is not None:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')


def get_arr(prev, curr):
    delta = tuple((cj - pj) / 5 for pj, cj in zip(prev,curr))
    orth = (delta[1], -delta[0], 0)
    if orth == (0,0,0):
        orth = (0,delta[2],-delta[1])
    l1 = tuple(cj - dj + oj for cj, dj, oj in zip(curr, delta, orth))
    l2 = tuple(cj - dj - oj for cj, dj, oj in zip(curr, delta, orth))
    return [l1,curr,l2,curr]


def draw_lines(ax, points, arrow=False, **kwargs):
    dim = len(points[0])
    if arrow:
        points = points + get_arr(points[-2], points[-1])
    args = [[x[j] for x in points] for j in range(dim)]
    ax.plot(*args, **kwargs)


def _get_abs_point(point, shift, scale):
    return tuple((pj + sj) * scale for pj, sj in zip(point, shift))


def _get_proto_center_line(proto):
    center_line = []
    for cube in proto:
        center = tuple((cj + 0.5) / proto.div for cj in cube)
        center_line.append(center)
    return center_line
