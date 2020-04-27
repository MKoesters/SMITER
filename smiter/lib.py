"""Core functionality."""
import csv
from io import TextIOWrapper
from tempfile import _TemporaryFileWrapper
from typing import Dict

from smiter.params.default_params import default_mzml_params, default_peak_properties

PROTON = 1.00727646677


def calc_mz(mass: float, charge: int) -> float:
    """Calculate m/z.

    Args:
        mass (float): mass
        charge (int): charge

    Returns:
        float: mass to charge value
    """
    mass = float(mass)
    charge = int(charge)
    calc_mz = (mass + (charge * PROTON)) / charge
    return calc_mz


def check_mzml_params(mzml_params: dict) -> dict:
    """Summary.

    Args:
        mzml_params (dict): dict with mzml parameters

    Returns:
        dict: mzml parameters with default values set

    Raises:
        Exception: Fails im required parameter without default is missing
    """
    for default_param, default_value in default_mzml_params.items():
        # param not set and default param required
        if (mzml_params.get(default_param, None) is None) and (default_value is None):
            raise Exception(f"mzml parameter {default_param} is required by not set!")
        elif mzml_params.get(default_param, None) is None:
            mzml_params[default_param] = default_value
    return mzml_params


def check_peak_properties(peak_properties: dict) -> dict:
    """Summary.

    Args:
        peak_properties (dict): Dict containing molecule as key and property dict as value

    Returns:
        dict: Dict containing molecule as key and property dict as value

    Raises:
        Exception: Fails im required parameter without default is missing
    """
    for mol, properties in peak_properties.items():
        for default_param, default_value in default_peak_properties.items():
            if (properties.get(default_param, None) is None) and (
                default_value is None
            ):
                raise Exception(
                    f"mzml parameter {default_param} is required by not set!"
                )
            elif properties.get(default_param, None) is None:
                properties[default_param] = default_value
    return peak_properties


def csv_to_peak_properties(csv_file: str) -> Dict:
    """Write csv file from peak properties.

    Args:
        csv_file (str): Csv file with peak properties

    Returns:
        Dict: peak properties
    """
    peak_properties = {}
    with open(csv_file) as fin:
        reader = csv.DictReader(fin)
        for line_dict in reader:
            cc = line_dict["chemical_formula"]
            peak_properties[cc] = {
                "trivial_name": line_dict["trivial_name"],
                "charge": line_dict.get("charge", 2),
                "scan_start_time": float(line_dict["scan_start_time"]),
                # currently only gaussian peaks from csv
                "peak_function": "gauss",
                "peak_params": {"sigma": line_dict.get("sigma", 2)},
                "peak_scaling_factor": float(line_dict["peak_scaling_factor"]),
                "peak_width": line_dict.get("peak_width", 30),
            }
    return peak_properties


def peak_properties_to_csv(peak_properties: Dict[str, dict], csv_file: str) -> str:
    """Wrirte peak properties dict to csv file.

    Args:
        peak_properties (Dict[str, dict]): Dict containing molecule as key and property dict as value
        csv_file (file): name of the csv file

    Returns:
        str: name of the csv file
    """
    if not isinstance(csv_file, TextIOWrapper):
        csv_file = open(csv_file, "w")
    csv_filename = csv_file.name
    fieldnames = [
        "chemical_formula",
        "trivial_name",
        "charge",
        "scan_start_time",
        "peak_width",
        "peak_scaling_factor",
        "peak_function",
        "peak_params",
    ]
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    writer.writeheader()
    for cc, attribs in peak_properties.items():
        line = {
            "chemical_formula": cc,
            "trivial_name": peak_properties[cc].get("trivial_name", ""),
            "charge": peak_properties[cc].get("charge", 2),
            "scan_start_time": peak_properties[cc]["scan_start_time"],
            "peak_function": peak_properties[cc]["peak_function"],
            "peak_params": ",".join(
                [
                    f"{key}={val}"
                    for key, val in peak_properties[cc]["peak_params"].items()
                ]
            ),
            "peak_scaling_factor": peak_properties[cc].get("peak_scaling_factor", 1e3),
            "peak_width": peak_properties[cc]["peak_width"],
        }
        writer.writerow(line)
    csv_file.close()
    return csv_filename
