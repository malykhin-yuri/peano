from peano.curves import Curve
from peano.base_maps import BaseMap


#
# 2D curves
#

def get_hilbert_curve():
    """
    Example of fractal curve due to D.Hilbert.

    This curve has minimal known L_1 ratio (9).
    """
    pattern = ('jiJ', ['ji','ij','ij','JI'])  # chain code + bases
    return Curve.parse_basis([pattern])


def get_peano_curve():
    """First example of space-filling curve due to G.Peano."""
    chain = 'jji' + 'JJi' + 'jj'
    bases = ['ij','Ij','ij', 'iJ','IJ','iJ', 'ij','Ij','ij']
    pattern = (chain, bases)
    return Curve.parse_basis([pattern])


def get_scepin_bauman_curve():
    """
    Minimal 3*3 Peano Curve by E.V. Shchepin and K.E. Bauman.

    Proceedings of the Steklov Institute of Mathematics, 2008, Vol. 263, pp. 236--256.
    """
    proto = (  # as in peano curve
        (0, 0), (0, 1), (0, 2),
        (1, 2), (1, 1), (1, 0),
        (2, 0), (2, 1), (2, 2),
    )
    base_maps = [
        BaseMap.id_map(dim=2),
        BaseMap.parse('(x,y)->(1-y,x)'),  # rot(90)
        BaseMap.parse('(x,y)->(y,x)'),

        BaseMap.parse('(x,y)->(x,1-y)'),
        BaseMap.parse('(x,y)->(1-y,1-x)'),
        BaseMap.parse('(x,y)->(y,1-x)'),  # rot(-90)

        BaseMap.id_map(dim=2),
        BaseMap.parse('(x,y)->(1-y,x)'),  # rot(90)
        BaseMap.parse('(x,y)->(y,x)'),
    ]
    return Curve(dim=2, div=3, patterns=[(proto, base_maps)])


def get_peano5_curve():
    """5-div analog of original Peano curve."""
    id_map = BaseMap.id_map(2)
    x_map = BaseMap.parse('(x,y)->(1-x,y)')
    y_map = BaseMap.parse('(x,y)->(x,1-y)')
    xy_map = BaseMap.parse('(x,y)->(1-x,1-y)')
    proto = [
        (0, 0), (0, 1), (0, 2), (0, 3), (0, 4),
        (1, 4), (1, 3), (1, 2), (1, 1), (1, 0),
        (2, 0), (2, 1), (2, 2), (2, 3), (2, 4),
        (3, 4), (3, 3), (3, 2), (3, 1), (3, 0),
        (4, 0), (4, 1), (4, 2), (4, 3), (4, 4),
    ]
    base_maps = [
        id_map, x_map, id_map, x_map, id_map,
        y_map, xy_map, y_map, xy_map, y_map,
        id_map, x_map, id_map, x_map, id_map,
        y_map, xy_map, y_map, xy_map, y_map,
        id_map, x_map, id_map, x_map, id_map,
    ]
    return Curve(dim=2, div=5, patterns=[(proto, base_maps)])


def get_meurthe_curve():
    """Meurthe curve, equivalent to Schepin-Bauman curve."""
    pattern = ('jjiJJijj', ['ji','jI','ij','Ji','JI','iJ','ji','jI','ij'])
    return Curve.parse_basis([pattern])

    
def get_coil_curve():
    """Coil 2D 3-div curve, see Haverkort & Walderveen."""
    pattern = ('jjiJJijj', ['ji','jI','ji','Ji','JI','Ji','ji','jI','ji'])
    return Curve.parse_basis([pattern])


def get_serpentine_curve():
    """Serpentine 2D 3-div curve, see Haverkort & Walderveen."""
    pattern = ('jjiJJijj', ['ij','jI','ji','iJ','JI','iJ','ji','jI','ij'])
    return Curve.parse_basis([pattern])


def get_r_curve():
    """R-curve, 2D 3-div, see Haverkort & Walderveen."""
    pattern = ('jjiiJIJi', ['ji','ji','ij','ij','ij','IJ','JI','JI','ij'])
    return Curve.parse_basis([pattern])


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
    omega_pattern = ('jiJ', ['1iJ','1jI','1ji~','1IJ~'])
    beta_pattern = ('jiJ', ['1iJ','1jI','1ji~','0Ji'])
    return Curve.parse_basis([omega_pattern, beta_pattern])


def get_ARW_Curve():
    """
    AR^2W^2 tetrafractal curve.

    See reference to Haverkort & Walderveen in get_beta_omega_curve doc
    """
    r_pattern = (['i','Ij','i'], ['3ij','1Ji~','2jI','1iJ'])  # pnum=0
    f_pattern = ('jiJ',          ['3ji','2Ij~','1ij','1JI'])  # pnum=1
    p_pattern = ('jiJ',          ['0ji','1jI','0Ji','1JI'])   # pnum=2
    g_pattern = ('jiJ',          ['0ij','2jI','0Ji','3jI~'])  # pnum=3
    return Curve.parse_basis([r_pattern, f_pattern, p_pattern, g_pattern])


#
# 3D curves
#

def get_haverkort_curve_a26():
    """
    3D Haverkort A26 curve, best known monofractal in l2.

    Monofractal curve with time reversal.
    "An inventory of three-dimensional Hilbert space-filling curves", Herman Haverkort
    https://arxiv.org/abs/1109.2323
    Curve A26.0010 1011.1011 0011, see p.10, p.15, p.18

    Properties:
    ratio: linf=12.4
    gate: (0,0,0)->(1,0,0)
    """
    pattern = ('jkJijKJ', ['Jki~', 'Kij~', 'kij', 'IKJ~', 'iKJ', 'kIj~', 'KIj', 'JIk'])
    return Curve.parse_basis([pattern])


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
    pattern = ('jkJijKJ', ['iKJ','jIK','jIk~','JkI','Jki~','jik','jiK~','kiJ~'])
    return Curve.parse_basis([pattern])


def get_tokarev_curve():
    """
    3D monofractal curve defined by Tokarev.

    Definition is taken from Haverkort's inverntory,
    Curve A26.0000 0000.0000 0000
    """
    p0 = ('jkJijKJ', ['jki','kij','kij','iJK','iJK','KIj','KIj','JkI'])
    return Curve.parse_basis([p0])


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
    p0 = ('jkJijKJ', ['1ijk','0jIK','1kJI','1JiK','1ijk','1jKI','1KJi','0JIk'])
    p1 = ('jkJiKjk', ['1ijk','0jIK','1kJI','1JiK','0ijk','1KjI','1jki','0kIJ'])
    return Curve.parse_basis([p0, p1])


def get_luna_curve():
    """
    3D bifractal Luna curve

    Inventory ref, p.16

    gates: (0,0,0)->(1,0,0), (0,0,0)->(1,1,1)
    """
    p0 = ('jkJijKJ', ['1kji','0jIK','1JIk','1iKJ','1ijk','1jKI','0JiK','1KiJ'])
    p1 = ('jkJiKjk', ['1kji','0jIK','1JIk','1iKJ','0ijk','1KjI','0kij','1jik'])
    return Curve.parse_basis([p0, p1])


def get_17_curve():
    """
    3D bifractal facet-gated curve with l2-ratio <17
    """
    p0 = ('jkiKJkI', ['1JKI~','0jKI','1kji','0kiJ~','1KiJ','0JKi','1kJI','1IkJ'])
    p1 = ('jiJkjIJ', ['1JKI~','0Ijk~','0jiK~','0KJI~','0Jki~','0ijK~','0IjK','1JIk'])
    return Curve.parse_basis([p0, p1])
