from ._nurbs import *
import numpy as np


def nurbs_base_from_freecad(bs):
    deg_u = bs.UDegree
    deg_v = bs.VDegree
    knots_u = np.array(bs.UKnotSequence)
    knots_v = np.array(bs.VKnotSequence)
    weights = np.array([j for i in bs.getWeights() for j in i])
    return NurbsBase2D(knots_u, knots_v, weights, deg_u, deg_v)


def uv_tesselation(bs, num_u=10, num_v=10):
    u1 = bs.getUKnot(bs.FirstUKnotIndex)
    u2 = bs.getUKnot(bs.LastUKnotIndex)
    v1 = bs.getVKnot(bs.FirstVKnotIndex)
    v2 = bs.getVKnot(bs.LastVKnotIndex)
    u_space = np.linspace(u1, u2, num_u)
    v_space = np.linspace(v1, v2, num_v)
    return np.array([[u, v] for u in u_space for v in v_space])


def pole_array(bs):
    poles = np.array([j for i in bs.getPoles() for j in i])
    return poles


def nurbs_base_1(bs):
    deg = bs.Degree
    knots = np.array(bs.KnotSequence)
    weights = np.array(bs.getWeights())
    return NurbsBase1D(knots, weights, deg)


def tesselation_1(bs, num=10):
    u1 = bs.getKnot(bs.FirstUKnotIndex)
    u2 = bs.getKnot(bs.LastUKnotIndex)
    return np.linspace(u1, u2, num)


def poles_1(bs):
    poles = np.array(bs.getPoles())
    return poles
