import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import binom
from typing import List


def model_C16(param, n):
    """_summary_

    Args:
        param: _description_
        n: _description_

    Returns:
        _description_
    """
    x = np.arange(0, n)
    uptake = param[0]
    p = param[1]
    density = binom.pmf(x, n-1, p)
    dist = np.zeros(n)
    dist[0] += uptake
    dist += (1 - uptake) * density
    return dist


def model_C18(param, n):
    """_summary_

    Args:
        param: _description_
        n: _description_

    Returns:
        _description_
    """
    x = np.arange(0, n)
    uptake = param[0]
    uptake_C16 = param[1]
    p = param[2]
    density = binom.pmf(x, n-1, p)
    dist = np.zeros(n)
    dist[0] += uptake
    dist[1] += (1 - uptake) * uptake_C16 
    dist += (1 - uptake) * (1 - uptake_C16) * density
    return dist


def square_dist(param: List, iso_even, model):
    """Calculate distance between measured distribution and prediction of the model with given parameters

    Args:
        param: vector of parameters of the model
        iso_even: vector with *only even* peaks of the measured isotopologue distribution
        model: function which predicts intensity of the even peaks in isotopologue distribution (model_C16 or model_C18 in this case)

    Returns:
        Sum of squared difference between experimental and predicted value for even peaks of isotopologue distribution
    """
    n = iso_even.shape[0]
    dist = model(param, n)
    return np.sum(np.square(iso_even - dist))


def fit_binomial(iso_dist: np.array, model: str):
    """Fit binomial model as described as Tumanov et al to the isotopologue distribution of C13-labeled fatty acids

    Args:
        iso_dist: isotopologue distribution - np array with normalized and corrected intensity of each isotopologe including M+0
        model: C16 or C18

    Returns:
        x - numbers of isotopologues
        n - number of even-numbered peaks
        param - parameters of the model
        success - True if minimize reached convergence condition
    """
    iso_even = iso_dist[0::2]
    n = iso_even.shape[0]
    x = np.arange(0, n)
    if model == "C16":
        param_min = minimize(square_dist, [iso_even[0], iso_even[1:].argmax() / n], args=(iso_even, model_C16), bounds=((0, 1), (0, 1)))
        # print(param_min)
        return x, n, param_min.x, param_min.success
    else:
        param_min = minimize(square_dist, [iso_even[0], iso_even[1] / iso_even[0], iso_even[2:].argmax() / n], args=(iso_even, model_C18), bounds=((0, 1), (0, 1), (0, 1)))
        return x, n, param_min.x, param_min.success

    
def calculate_mean(iso_dist_even, model_type="C16"):
    """Calculate the position of the mean of the isotopologue distribution

    Args:
        iso_dist_even: vector with *only even* peaks of the measured isotopologue distribution

    Returns:
        mean - mean of the isotopologue distribution (2 * I(M+2) + 4 * I(M+4)... / sum(I))
        p - parameter p of the binomial distribution estimated as mean / n
    """
    iso_dist_model = iso_dist_even.copy()
    iso_dist_model[0] = 0
    if model_type == "C18":
        iso_dist_model[1] = 0
        
    n = len(iso_dist_model) - 1
    x = np.arange(0, len(iso_dist_model))
    iso_dist_model = (iso_dist_model / (np.sum(iso_dist_model) + 1e-7))
    mean_M = np.sum(x * iso_dist_model)
    p = mean_M / n
    return mean_M * 2, p


def get_full_prediction(dist: np.array, prediction: np.array):
    full_prediction = np.zeros_like(dist)
    full_prediction[::2] = prediction
    return full_prediction