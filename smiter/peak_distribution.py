#!/usr/bin/env python3
"""Distribution funtion for chromo peaks.

Attributes:
    distributions (dict): mapping distribution name to distribution function
"""
import math
from typing import Callable, Dict

from scipy.stats import gamma


def gauss_dist(rt: float, sigma: float = 1, mu: float = 0):
    """Calc Gauss distribution.

    Args:
        rt (float): retention time
        sigma (float, optional): standard deviation
        mu (float, optional): mean

    Returns:
        float: intensity from distribution at time `rt`
    """
    return (
        1
        / (sigma * math.sqrt(2 * math.pi))
        * pow(math.e, (-0.5 * pow(((rt - mu) / sigma), 2)))
    )


def gamma_dist(rt: float, a: float = 5, scale: float = 0.33):
    """Calc gamma distribution.

    Args:
        rt (float): retention time
        a (float, optional): Description
        scale (float, optional): Description

    Returns:
        float: intensity from distribution at time `rt`
    """
    return gamma.pdf(rt, a=a, scale=scale)


distributions: Dict[str, Callable] = {
    "gauss": gauss_dist,
    "gamma": gamma_dist,
}
