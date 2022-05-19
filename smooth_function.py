import matplotlib.pyplot as plt
from math import *
import numpy as np

def smooth(step_curve, j_max=None, width_max=None) :
    if j_max == None :
        j_max = max(0,int(len(step_curve)/100))
    if width_max == None :
        width_max = (np.max(step_curve)-np.min(step_curve))/2
    curve_copy = np.copy(step_curve)
    for i in range(len(step_curve)) :
        count = 0
        sum = 0
        var_j = j_max
        while not((0 <= i+var_j < len(step_curve)) and (0 <= i-var_j < len(step_curve))) :
            var_j -= 1
        while not((np.abs(step_curve[i+var_j]-step_curve[i]) < width_max) and (np.abs(step_curve[i-var_j]-step_curve[i]) < width_max)) :
            var_j -= 1
        for j in range(-var_j,var_j+1) :
            count += 1+var_j-np.abs(j)
            sum += step_curve[i+j]*(1+var_j-np.abs(j))
        if var_j != 0 :
            curve_copy[i] = sum/count
    return curve_copy

def correct(smoothed_curve, init_curve, j_max=None) :
    if j_max == None :
        j_max = max(0,int(len(smoothed_curve)/100))
    list_mult = np.empty(len(smoothed_curve))
    for i in range(len(smoothed_curve)) :
        var_j = j_max
        while not((0 <= i+var_j < len(smoothed_curve)) and (0 <= i-var_j < len(smoothed_curve))) :
            var_j -= 1
        mult = np.sum(smoothed_curve)/len(smoothed_curve)
        fact = mult/10
        precision = mult/10000
        above = (np.sum(smoothed_curve[i-var_j:i+var_j+1]) > np.sum(init_curve[i-var_j:i+var_j+1]))
        while np.abs(np.sum(smoothed_curve[i-var_j:i+var_j+1])*mult-np.sum(init_curve[i-var_j:i+var_j+1])) > precision :
            if np.sum(smoothed_curve[i-var_j:i+var_j+1])*mult > np.sum(init_curve[i-var_j:i+var_j+1]) :
                if not above :
                    fact = fact*0.1
                    above = True
                mult -= fact
            else :
                if above :
                    fact = fact*0.1
                    above = False
                mult += fact
        list_mult[i] = mult
    return smooth(list_mult)

def recur_smooth(step_curve, init_curve, nb_recur=5) :
    if nb_recur == 0 :
        return step_curve
    else :
        smoothed_curve = smooth(step_curve)
        correction = correct(smoothed_curve, init_curve, j_max=int(len(smoothed_curve)/100+10*(1/5*(1-nb_recur)+1)))
        return recur_smooth(smoothed_curve*correction, init_curve, nb_recur-1)

def recur_smooth2(step_curve, init_curve, nb_recur=5) :
    if nb_recur == 0 :
        return step_curve
    else :
        smoothed_curve = smooth(step_curve)
        return recur_smooth2(smoothed_curve, init_curve, nb_recur-1)

def final_smooth(step_curve, norm=False) :
    really_smoothed_curve = recur_smooth(step_curve, smooth(step_curve))
    if norm :
        mult = np.sum(smoothed_curve)/len(smoothed_curve)
        fact = mult/10
        precision = mult/10000
        above = (np.sum(smoothed_curve) > 1)
        while np.abs(np.sum(smoothed_curve)*mult-1) > precision :
            if np.sum(smoothed_curve)*mult > 1 :
                if not above :
                    fact = fact*0.1
                    above = True
                mult -= fact
            else :
                if above :
                    fact = fact*0.1
                    above = False
                mult += fact
        return really_smoothed_curve*mult
    else :
        return really_smoothed_curve

def final_smooth2(step_curve, norm=False) :
    really_smoothed_curve = recur_smooth2(step_curve, smooth(step_curve))
    if norm :
        mult = np.sum(smoothed_curve)/len(smoothed_curve)
        fact = mult/10
        precision = mult/10000
        above = (np.sum(smoothed_curve) > 1)
        while np.abs(np.sum(smoothed_curve)*mult-1) > precision :
            if np.sum(smoothed_curve)*mult > 1 :
                if not above :
                    fact = fact*0.1
                    above = True
                mult -= fact
            else :
                if above :
                    fact = fact*0.1
                    above = False
                mult += fact
        return really_smoothed_curve*mult
    else :
        return really_smoothed_curve
