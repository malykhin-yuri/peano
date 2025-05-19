"""
Known Peano curves with some information about them.
"""

from dataclasses import dataclass
from typing import Any

from quicktions import Fraction  # type: ignore

from peano.curves import Curve
from peano.subsets import Link, Point


# TODO: better typing for DilationInfo (use enum MT & dict MT->Info?; spec - ?)
@dataclass
class DilationInfo:
    type DilationValue = Fraction | tuple[float, float] | None
    l1: DilationValue = None
    l2: DilationValue = None
    l2_squared: DilationValue = None
    linf: DilationValue = None
    spec: dict[str, Any] | None = None


@dataclass
class CurveInfo:
    curve: Curve
    dilation: DilationInfo | None = None
    moments: dict[int, list[Fraction]] | None = None
    gates: list[Link[Point]] | None = None
    depth: int | None = None
    junctions_count: int | None = None


def get_all_curves() -> list[CurveInfo]:
    return [
        get_hilbert_curve(),
        get_peano_curve(),
        get_bauman_curve(),
        get_scepin_bauman_curve(),
        get_peano5_curve(),
        get_meurthe_curve(),
        get_coil_curve(),
        get_serpentine_curve(),
        get_r_curve(),
        get_ye_curve(),
        get_beta_omega_curve(),
        get_ARW_Curve(),
        get_haverkort_curve_a26(),
        get_haverkort_curve_f(),
        get_tokarev_curve(),
        get_neptunus_curve(),
        get_luna_curve(),
        get_iupiter_curve(),
        get_spring_curve(),
        get_3d_l1_best_curve(),
        get_4d_facet_gated_curve(),
    ]


#
# 2D curves
#

def get_hilbert_curve() -> CurveInfo:
    """
    Example of fractal curve due to D.Hilbert.
    This curve has minimal known L_1 ratio.
    """
    pattern = ('jiJ', 'ji,ij,ij,JI')
    return CurveInfo(
        curve=Curve.parse([pattern]),
        dilation=DilationInfo(l1=9, l2=6, linf=6),
        depth=2,
    )


def get_peano_curve() -> CurveInfo:
    """First example of space-filling curve due to G.Peano."""
    pattern = ('jjiJJijj', 'ij,Ij,ij,iJ,IJ,iJ,ij,Ij,ij')
    return CurveInfo(
        curve=Curve.parse([pattern]),
        dilation=DilationInfo(l1=Fraction(32, 3), l2=8, linf=8),
        depth=1,
    )


def get_bauman_curve() -> CurveInfo:
    """
    TODO: add description
    """
    pattern = ('ijIjiiJJ', 'ij,ji,Ji~,iJ~,ij,ij,Ij~,jI~,jI~')
    return CurveInfo(
        curve=Curve.parse([pattern]),
        dilation=DilationInfo(l2=Fraction(17, 3), linf=Fraction(9, 2)),
    )


def get_scepin_bauman_curve() -> CurveInfo:
    """
    Minimal 3*3 Peano Curve by E.V. Shchepin and K.E. Bauman.

    Proceedings of the Steklov Institute of Mathematics, 2008, Vol. 263, pp. 236--256.
    """
    pattern = ('jjiJJijj', 'ij,jI,ji,iJ,JI,Ji,ij,jI,ji')
    return CurveInfo(
        curve=Curve.parse([pattern]),
        dilation=DilationInfo(l1=Fraction(32, 3), l2=Fraction(17, 3), linf=Fraction(16, 3)),
    )


def get_peano5_curve() -> CurveInfo:
    """5-div analog of original Peano curve."""
    pattern = (
        'jjjjiJJJJijjjjiJJJJijjjj',
        'ij,Ij,ij,Ij,ij, iJ,IJ,iJ,IJ,iJ, ij,Ij,ij,Ij,ij, iJ,IJ,iJ,IJ,iJ, ij,Ij,ij,Ij,ij',
    )
    return CurveInfo(
        curve=Curve.parse([pattern]),
    )


def get_meurthe_curve() -> CurveInfo:
    """Meurthe curve, equivalent to Schepin-Bauman curve."""
    pattern = ('jjiJJijj', 'ji,jI,ij,Ji,JI,iJ,ji,jI,ij')
    return CurveInfo(
        curve=Curve.parse([pattern]),
        dilation=DilationInfo(l1=Fraction(32, 3), l2=Fraction(17, 3), linf=Fraction(16, 3)),
    )

    
def get_coil_curve() -> CurveInfo:
    """Coil 2D 3-div curve, see Haverkort & Walderveen."""
    pattern = ('jjiJJijj', 'ji,jI,ji,Ji,JI,Ji,ji,jI,ji')
    return CurveInfo(
        curve=Curve.parse([pattern]),
        dilation=DilationInfo(l1=Fraction(32, 3), l2=Fraction(20, 3), linf=Fraction(20, 3)),
    )


def get_serpentine_curve() -> CurveInfo:
    """Serpentine 2D 3-div curve, see Haverkort & Walderveen."""
    pattern = ('jjiJJijj', 'ij,jI,ji,iJ,JI,iJ,ji,jI,ij')
    return CurveInfo(
        curve=Curve.parse([pattern]),
        dilation=DilationInfo(l1=10, l2=Fraction(25, 4), linf=Fraction(45, 8)),
    )


def get_r_curve() -> CurveInfo:
    """R-curve, 2D 3-div, see Haverkort & Walderveen."""
    pattern = ('jjiiJIJi', 'ji,ji,ij,ij,ij,IJ,JI,JI,ij')
    return CurveInfo(
        curve=Curve.parse([pattern]),
        dilation=DilationInfo(l1=Fraction(32, 3), l2=Fraction(20, 3), linf=Fraction(20, 3)),
    )


def get_ye_curve() -> CurveInfo:
    """YE-curve: 5*5 monofractal with l2-ratio 5 43/73."""
    pattern = (
        'jiJijjIIjjiJijiiJIJiJIJi',
        '0Ji~,0ij,0ij,0jI~,0Ji~,0Ji~,0Ji~,0IJ,0IJ,0Ji~,0ij,0ij,0jI~,0Ji~,0ij,0ij,0ij,0IJ,0jI~,0jI~,0ij,0IJ,0jI~,0jI~,0ij',
    )
    return CurveInfo(
        curve=Curve.parse([pattern]),
        dilation=DilationInfo(l2=Fraction(408, 73)),
        depth=1,
    )


# 2D polyfractals

def get_beta_omega_curve() -> CurveInfo:
    """
    Beta-Omega bifractal 2d curve, best polyfractal in l2

    See J.-M. Wierum, ``Definition of a New Circular Space-Filling Curve'',
    Technical report PC^2, 2002

    Realization from Haverkort & Walderveen,
    ``Locality and bounding-box quality of two-dimensional space-filling curves'',
    Comp.Geom and Appl, v.43 (2010)
    """
    omega_pattern = ('jiJ', '1iJ,1jI,1ji~,1IJ~')
    beta_pattern = ('jiJ', '1iJ,1jI,1ji~,0Ji')
    return CurveInfo(
        curve=Curve.parse([omega_pattern, beta_pattern]),
        dilation=DilationInfo(l1=9, l2=5, linf=5),
        gates=[Link.parse_gates('(0,1/3)->(1,1/3)'), Link.parse_gates('(0,1/3)->(2/3,0)')],
    )


def get_ARW_Curve() -> CurveInfo:
    """
    AR^2W^2 tetrafractal curve.

    See reference to Haverkort & Walderveen in get_beta_omega_curve doc
    """
    r_pattern = ('i(Ij)i', '3ij,1Ji~,2jI,1iJ')  # pnum=0
    f_pattern = ('jiJ',    '3ji,2Ij~,1ij,1JI')  # pnum=1
    p_pattern = ('jiJ',    '0ji,1jI,0Ji,1JI')   # pnum=2
    g_pattern = ('jiJ',    '0ij,2jI,0Ji,3jI~')  # pnum=3
    return CurveInfo(
        curve=Curve.parse([r_pattern, f_pattern, p_pattern, g_pattern]),
        dilation=DilationInfo(l1=12, l2=Fraction(260, 43), linf=Fraction(27, 5)),  # l2,l2inf: Haverkort & Walderveen
    )


#
# 3D curves
#

def get_haverkort_curve_a26() -> CurveInfo:
    """
    3D Haverkort A26 curve, best monofractal in linf.

    Monofractal curve with time reversal.
    "An inventory of three-dimensional Hilbert space-filling curves", Herman Haverkort
    https://arxiv.org/abs/1109.2323
    Curve A26.0010 1011.1011 0011, see p.10, p.15, p.18
    """
    pattern = ('jkJijKJ', 'Jki~,Kij~,kij,IKJ~,iKJ,kIj~,KIj,JIk')
    return CurveInfo(
        curve=Curve.parse([pattern]),
        dilation=DilationInfo(l1=Fraction(99, 1) + Fraction(5, 9), l2=[22.7,22.9], linf=Fraction(12, 1) + Fraction(4, 9)),
        moments={0: [Fraction(k, 28) for k in [0, 5, 9, 12, 16, 19, 23, 28]]},
        gates=[Link.parse_gates('(0,0,0)->(1,0,0)')],
    )


def get_haverkort_curve_f() -> CurveInfo:
    """
    3D Haverkort F curve, best known monofractal in l2

    Monofractal curve with time reversal,
    "An inventory of three-dimensional Hilbert space-filling curves", Herman Haverkort
    https://arxiv.org/abs/1109.2323
    Curve F, see p.13, p.15, p.18

    Properties:
    l2 ratio is best known among monofractals!
    """
    pattern = ('jkJijKJ', 'iKJ,jIK,jIk~,JkI,Jki~,jik,jiK~,kiJ~')
    return CurveInfo(
        curve=Curve.parse([pattern]),
        dilation=DilationInfo(l1=[89.754, 89.758], l2=[18.5,18.7], linf=14),
        moments={0: [Fraction(k, 28) for k in [1, 6, 8, 13, 15, 20, 22, 27]]},
        gates=[Link.parse_gates('(0,1/3,1/3)->(2/3,1/3,0)')],
    )


def get_tokarev_curve() -> CurveInfo:
    """
    3D monofractal curve defined by Tokarev.

    Definition is taken from Haverkort's inventory,
    Curve A26.0000 0000.0000 0000 (page 9, Fig.5(b))
    """
    p0 = ('jkJijKJ', 'jki,kij,kij,iJK,iJK,KIj,KIj,JkI')
    return CurveInfo(
        curve=Curve.parse([p0]),
        dilation=DilationInfo(
            l1=[98.2, 98.4],
            l2_squared=Fraction(5215408884, 7579009),  # TODO ref
            linf=Fraction(896, 37),
            spec={'linf': {'face_dim': 2}},
        ),
        # from Scepin & Korneev (2018):
        moments={
            0: [Fraction(k, 126) for k in [0, 22, 41, 50, 76, 85, 104, 126]],
            1: [Fraction(k, 4194176) for k in [0, 693632, 1364617, 1659520]]
                + [Fraction(k, 65534) for k in [0, 11433, 38292, 44200]]
                + [Fraction(k, 524272) for k in [0, 169360, 316073, 431496]],
            -1: [Fraction(k, 4194176) for k in [2534656, 2829559, 3500544, 4194176]]
                + [Fraction(k, 65534) for k in [21334, 27242, 54101, 65534]]
                + [Fraction(k, 524272) for k in [92776, 208199, 354912, 524272]],
            2: [Fraction(k, 224) for k in [0, 0, 0, 37, 72, 128]],
            -2: [Fraction(k, 224) for k in [96, 152, 187, 224, 224, 224]],
        },
        depth=3,
        junctions_count=9,
    )


# 3D polyfractals

def get_neptunus_curve() -> CurveInfo:
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
    return CurveInfo(
        curve=Curve.parse([p0, p1]),
        dilation=DilationInfo(l1=[88.8, 89.0], l2=[18.2, 18.4], linf=Fraction(945, 100)),
        gates=[Link.parse_gates('(0,0,0)->(1,0,0)'), Link.parse_gates('(0,0,0)->(1,1,1)')],
    )


def get_luna_curve() -> CurveInfo:
    """
    3D bifractal Luna curve

    Inventory ref, p.16
    """
    p0 = ('jkJijKJ', '1kji,0jIK,1JIk,1iKJ,1ijk,1jKI,0JiK,1KiJ')
    p1 = ('jkJiKjk', '1kji,0jIK,1JIk,1iKJ,0ijk,1KjI,0kij,1jik')
    return CurveInfo(
        curve=Curve.parse([p0, p1]),
        dilation=DilationInfo(l1=[75.5, 75.7], l2=[18.2, 18.4], linf=14),
        gates=[Link.parse_gates('(0,0,0)->(1,0,0)'), Link.parse_gates('(0,0,0)->(1,1,1)')],
    )


def get_iupiter_curve() -> CurveInfo:
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

    return CurveInfo(
        curve=Curve.parse([pa0, pb1, pc2, pd3, pe4]),
        dilation=DilationInfo(l1=[88.6, 88.8], l2=[24.8, 30.0], linf=[16.9, 17.1]),
        gates=[
            Link.parse_gates('(0,2/5,1/5)->(4/5,2/5,0)'),
            Link.parse_gates('(0,2/5,1/5)->(4/5,0,2/5)'),
            Link.parse_gates('(0,2/5,1/5)->(1,2/5,1/5)'),
            Link.parse_gates('(0,1/5,2/5)->(4/5,2/5,0)'),
            Link.parse_gates('(0,2/5,1/5)->(4/5,1,3/5)'),
        ],
    )



def get_spring_curve() -> CurveInfo:
    """Spring curve - 3D bifractal facet-gated curve with l2-ratio <17"""
    # this curve coincides with the spring example in the Scepin-Malykhin paper
    # calculations show that WD = 1533\sqrt{6}/221
    p0 = ('jikIJiK', '1KIJ~,0KIj,1kji,0Jki~,1JkI,0kIJ,1KJi,1JiK')
    p1 = ('jkJijKJ', '1KIJ~,0ijK~,0Ikj~,0KJI~,0kiJ~,0Ijk~,0IjK,1iKJ')
    return CurveInfo(
        curve=Curve.parse([p0, p1]),
        dilation=DilationInfo(l1=[82.9, 83.0], l2=[16.9, 17.0]),
    )


def get_3d_l1_best_curve() -> CurveInfo:
    """
    3d-curve with best l1-dilation 89.74, found by:
    $ python search.py --dim 3 --mult 1 --div 2 --gates '(0,0,0)->(0,0,1)' --metric l1 --rel-tol-inv 100000
    """
    p0 = ('jiJkjIJ', 'ikJ~,kjI~,kjI~,JIk,JIK~,KjI,jKi~,iKj~')
    return CurveInfo(
        curve=Curve.parse([p0]),
        dilation=DilationInfo(l1=[89.742, 89.745]),
    )


def get_4d_facet_gated_curve() -> CurveInfo:
    """4d facet-gated curve with dilation < 61.9354839"""
    # WD = 1920/31?
    p = ('kljKLkiKlkJLKlI', 'iKLJ,kLij,lJki,jklI,lKjI~,LKIj,IkLj~,ikLj,LKij~,lKji,kjlI,lJkI~,LkIJ,LKji~,iljk~,jIlK~')
    return CurveInfo(
        curve=Curve.parse([p]),
    )
