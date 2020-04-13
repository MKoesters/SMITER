"""Callables for injection noise into scans.

Upon calling the callabe, a list/np.array of mz and intensities should be returned.
Arguments should be passed via *args and **kwargs
"""
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple

import numpy as np
import pyqms

from smiter.lib import calc_mz
from smiter.synthetic_mzml import Scan

easter_egg = """

                   .-^-.
                 .'=^=^='.
                /=^=^=^=^=\

        .-~-.  :^= HAPPY =^;
      .'~~*~~'.|^ EASTER! ^|
     /~~*~~~*~~\^=^=^=^=^=^:
    :~*~~~*~~~*~;\.-*))`*-,/
    |~~~*~~~*~~|/*  ((*   *'.
    :~*~~~*~~~*|   *))  *   *\

     \~~*~~~*~~| *  ((*   *  /
      `.~~*~~.' \  *))  *  .'
        `~~~`    '-.((*_.-'
"""


class AbstractNoiseInjector(ABC):
    """Summary."""

    def __init__(self, *args, **kwargs):
        """Initialize noise injector."""
        pass  # pragma: no cover

    @abstractmethod
    def inject_noise(self, scan: Scan, *args, **kwargs):
        """Main noise injection method.

        Args:
            scan (Scan): Scan object
            *args: Description
            **kwargs: Description

        """
        pass  # pragma: no cover

    @abstractmethod
    def _ms1_noise(self, scan: Scan, *args, **kwargs):
        pass

    @abstractmethod
    def _msn_noise(self, scan: Scan, *args, **kwargs):
        pass


class TestNoiseInjector(AbstractNoiseInjector):
    def __init__(self, *args, **kwargs):
        print(easter_egg)
        self.args = args
        self.kwargs = kwargs

    def inject_noise(self, scan: Scan, *args, **kwargs):
        """Main noise injection method.

        Args:
            scan (Scan): Scan object
            *args: Description
            **kwargs: Description

        """
        if scan.ms_level == 1:
            scan = self._ms1_noise(scan)
        elif scan.ms_level > 1:
            scan = self._msn_noise(scan)
        return scan

    def _ms1_noise(self, scan: Scan, *args, **kwargs):
        """Generate ms1 noise.

        Args:
            scan (Scan): Description
            *args: Description
            **kwargs: Description
        """
        mz_noise = self._generate_mz_noise(scan, *args, **kwargs)
        intensity_noise = self._generate_intensity_noise(scan, *args, **kwargs)
        scan.mz += mz_noise
        scan.i += intensity_noise
        return scan

    def _msn_noise(self, scan: Scan, *args, **kwargs):
        """Generate msn noise.

        Args:
            scan (Scan): Description
            *args: Description
            **kwargs: Description
        """
        mz_noise = self._generate_mz_noise(scan, *args, **kwargs)
        intensity_noise = self._generate_intensity_noise(scan, *args, **kwargs)
        scan.mz += mz_noise
        scan.i += intensity_noise
        return scan

    def _generate_mz_noise(self, scan: Scan, *args, **kwargs):
        """Generate noise for mz_array.

        Args:
            scan (Scan): Scan object
            *args: Description
            **kwargs: Description
        """
        noise = [0 for i in range(len(scan.mz))]
        return noise

    def _generate_intensity_noise(
        self, scan: Scan, *args, **kwargs
    ):
        """Generate intensity noise.

        Args:
            scan (Scan): Scan object
            *args: Description
            **kwargs: Description
        """
        np.random.seed(1312)
        noise = np.random.normal(
            max(scan.i) * kwargs.get("max_noise_level", -0.05),
            max(scan.i) * kwargs.get("max_noise_level", 0.05),
            len(scan.i)
        )
        # breakpoint()
        return noise
