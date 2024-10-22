from peano.curves import Curve


#
# 2D curves
#

def get_hilbert_curve():
    """
    Example of fractal curve due to D.Hilbert.

    This curve has minimal known L_1 ratio (9).
    """
    pattern = ('jiJ', 'ji,ij,ij,JI')
    return Curve.parse([pattern])


def get_peano_curve():
    """First example of space-filling curve due to G.Peano."""
    pattern = ('jjiJJijj', 'ij,Ij,ij,iJ,IJ,iJ,ij,Ij,ij')
    return Curve.parse([pattern])


def get_bauman_curve():
    """
    TODO: add description
    """
    pattern = ('ijIjiiJJ', 'ij,ji,Ji~,iJ~,ij,ij,Ij~,jI~,jI~')
    return Curve.parse([pattern])


def get_scepin_bauman_curve():
    """
    Minimal 3*3 Peano Curve by E.V. Shchepin and K.E. Bauman.

    Proceedings of the Steklov Institute of Mathematics, 2008, Vol. 263, pp. 236--256.
    """
    pattern = ('jjiJJijj', 'ij,jI,ji,iJ,JI,Ji,ij,jI,ji')
    return Curve.parse([pattern])


def get_peano5_curve():
    """5-div analog of original Peano curve."""
    pattern = (
        'jjjjiJJJJijjjjiJJJJijjjj',
        'ij,Ij,ij,Ij,ij, iJ,IJ,iJ,IJ,iJ, ij,Ij,ij,Ij,ij, iJ,IJ,iJ,IJ,iJ, ij,Ij,ij,Ij,ij',
    )
    return Curve.parse([pattern])


def get_meurthe_curve():
    """Meurthe curve, equivalent to Schepin-Bauman curve."""
    pattern = ('jjiJJijj', 'ji,jI,ij,Ji,JI,iJ,ji,jI,ij')
    return Curve.parse([pattern])

    
def get_coil_curve():
    """Coil 2D 3-div curve, see Haverkort & Walderveen."""
    pattern = ('jjiJJijj', 'ji,jI,ji,Ji,JI,Ji,ji,jI,ji')
    return Curve.parse([pattern])


def get_serpentine_curve():
    """Serpentine 2D 3-div curve, see Haverkort & Walderveen."""
    pattern = ('jjiJJijj', 'ij,jI,ji,iJ,JI,iJ,ji,jI,ij')
    return Curve.parse([pattern])


def get_r_curve():
    """R-curve, 2D 3-div, see Haverkort & Walderveen."""
    pattern = ('jjiiJIJi', 'ji,ji,ij,ij,ij,IJ,JI,JI,ij')
    return Curve.parse([pattern])


def get_ye_curve():
    """YE-curve: 5*5 monofractal with l2-ratio 5 43/73."""
    pattern = (
        'jiJijjIIjjiJijiiJIJiJIJi',
        '0Ji~,0ij,0ij,0jI~,0Ji~,0Ji~,0Ji~,0IJ,0IJ,0Ji~,0ij,0ij,0jI~,0Ji~,0ij,0ij,0ij,0IJ,0jI~,0jI~,0ij,0IJ,0jI~,0jI~,0ij',
    )
    return Curve.parse([pattern])


# 2D polyfractals

def get_beta_omega_curve():
    """
    Beta-Omega bifractal 2d curve, best polyfractal in l2

    See J.-M. Wierum, ``Definition of a New Circular Space-Filling Curve'',
    Technical report PC^2, 2002

    Realization from Haverkort & Walderveen,
    ``Locality and bounding-box quality of two-dimensional space-filling curves'',
    Comp.Geom and Appl, v.43 (2010)

    ratio: l2=5, l1=9, linf=5
    gates: omega: (0,1/3)->(1,1/3), beta: (0,1/3)->(2/3,0)
    """
    omega_pattern = ('jiJ', '1iJ,1jI,1ji~,1IJ~')
    beta_pattern = ('jiJ', '1iJ,1jI,1ji~,0Ji')
    return Curve.parse([omega_pattern, beta_pattern])


def get_ARW_Curve():
    """
    AR^2W^2 tetrafractal curve.

    See reference to Haverkort & Walderveen in get_beta_omega_curve doc
    """
    r_pattern = ('i(Ij)i', '3ij,1Ji~,2jI,1iJ')  # pnum=0
    f_pattern = ('jiJ',    '3ji,2Ij~,1ij,1JI')  # pnum=1
    p_pattern = ('jiJ',    '0ji,1jI,0Ji,1JI')   # pnum=2
    g_pattern = ('jiJ',    '0ij,2jI,0Ji,3jI~')  # pnum=3
    return Curve.parse([r_pattern, f_pattern, p_pattern, g_pattern])


#
# 3D curves
#

def get_haverkort_curve_a26():
    """
    3D Haverkort A26 curve, best monofractal in linf.

    Monofractal curve with time reversal.
    "An inventory of three-dimensional Hilbert space-filling curves", Herman Haverkort
    https://arxiv.org/abs/1109.2323
    Curve A26.0010 1011.1011 0011, see p.10, p.15, p.18

    Properties:
    ratio: linf=12.4
    gate: (0,0,0)->(1,0,0)
    """
    pattern = ('jkJijKJ', 'Jki~,Kij~,kij,IKJ~,iKJ,kIj~,KIj,JIk')
    return Curve.parse([pattern])


def get_haverkort_curve_f():
    """
    3D Haverkort F curve, best known monofractal in l2

    Monofractal curve with time reversal,
    "An inventory of three-dimensional Hilbert space-filling curves", Herman Haverkort
    https://arxiv.org/abs/1109.2323
    Curve F, see p.13, p.15, p.18

    Properties:
    ratio: l1=89.8, l2=18.6 (best known!)
    gate: (0,1/3,1/3)->(2/3,1/3,0)
    """
    pattern = ('jkJijKJ', 'iKJ,jIK,jIk~,JkI,Jki~,jik,jiK~,kiJ~')
    return Curve.parse([pattern])


def get_tokarev_curve():
    """
    3D monofractal curve defined by Tokarev.

    Definition is taken from Haverkort's inventory,
    Curve A26.0000 0000.0000 0000 (page 9, Fig.5(b))
    """
    p0 = ('jkJijKJ', 'jki,kij,kij,iJK,iJK,KIj,KIj,JkI')
    return Curve.parse([p0])


# 3D polyfractals

def get_neptunus_curve():
    """
    3D bifractal Neptunus curve, best in linf (9.45)

    "An inventory of three-dimensional Hilbert space-filling curves", Herman Haverkort
    https://arxiv.org/abs/1109.2323
    p.16

    Properties:
    ratio: linf=9.45, best in class of poly-Hilbert curves
    gates: (0,0,0)->(1,0,0), (0,0,0)->(1,1,1)
    """
    p0 = ('jkJijKJ', '1ijk,0jIK,1kJI,1JiK,1ijk,1jKI,1KJi,0JIk')
    p1 = ('jkJiKjk', '1ijk,0jIK,1kJI,1JiK,0ijk,1KjI,1jki,0kIJ')
    return Curve.parse([p0, p1])


def get_luna_curve():
    """
    3D bifractal Luna curve

    Inventory ref, p.16

    gates: (0,0,0)->(1,0,0), (0,0,0)->(1,1,1)
    """
    p0 = ('jkJijKJ', '1kji,0jIK,1JIk,1iKJ,1ijk,1jKI,0JiK,1KiJ')
    p1 = ('jkJiKjk', '1kji,0jIK,1JIk,1iKJ,0ijk,1KjI,0kij,1jik')
    return Curve.parse([p0, p1])


def get_iupiter_curve():
    """Inventory, p.17"""
    pa0 = ('jkJijKJ', ['0ikJ', '4jkI', '1jki~', '0IkJ~', '0ikJ', '1jkI', '0jIK~', '0kIJ~'])

    pb1 = (pa0[0], pa0[1][:])
    pb1[1][-1] = '2JIK'

    pc2 = (pa0[0], pa0[1][:])
    pc2[1][-2] = '4jki~'
    pc2[1][-1] = '0IkJ~'

    pd3 = (pa0[0], pa0[1][:])
    pd3[1][0] = '4ijK'

    pe4 = ('jkJiKjk', pa0[1][:])
    pe4[1][-4] = '3iJk'
    pe4[1][-3] = '1KJI'
    pe4[1][-2] = '0KIj~'
    pe4[1][-1] = '0JIk~'

    return Curve.parse([pa0, pb1, pc2, pd3, pe4])


def get_spring_curve():
    """Spring curve - 3D bifractal facet-gated curve with l2-ratio <17"""
    # this curve coincides with the spring example in the Scepin-Malykhin paper
    # calculations show that WD = 1533\sqrt{6}/221
    p0 = ('jikIJiK', '1KIJ~,0KIj,1kji,0Jki~,1JkI,0kIJ,1KJi,1JiK')
    p1 = ('jkJijKJ', '1KIJ~,0ijK~,0Ikj~,0KJI~,0kiJ~,0Ijk~,0IjK,1iKJ')
    return Curve.parse([p0, p1])


def get_3d_l1_best_curve():
    """
    3d-curve with best l1-dilation 89.74, found by:
    $ python search.py --dim 3 --mult 1 --div 2 --gates '(0,0,0)->(0,0,1)' --metric l1 --rel-tol-inv 100000
    """
    p0 = ('jiJkjIJ', 'ikJ~,kjI~,kjI~,JIk,JIK~,KjI,jKi~,iKj~')
    return Curve.parse([p0])


def get_4d_facet_gated_curve():
    """4d facet-gated curve with dilation < 61.9354839"""
    # WD = 1920/31?
    p = ('kljKLkiKlkJLKlI', 'iKLJ,kLij,lJki,jklI,lKjI~,LKIj,IkLj~,ikLj,LKij~,lKji,kjlI,lJkI~,LkIJ,LKji~,iljk~,jIlK~')
    return Curve.parse([p])
