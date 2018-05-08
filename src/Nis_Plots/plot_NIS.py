#!/usr/bin/python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv

def plot_csv():

    with open("lidarNis.csv",'r')as fileReader:
        lidar_Nis = list(csv.reader(fileReader,delimiter=","))

    lidar_Nis_float =  list(map(float,lidar_Nis[0]))
    print("num elements:", len(lidar_Nis_float))
    lidar_95pct = [5.991]*len(lidar_Nis_float); 
    fig1 = plt.figure();
    ax1 = fig1.add_subplot(121);
    ax1.plot(lidar_Nis_float)
    ax1.plot(lidar_95pct)
    ax1.set_title("Lidar NIS")
    
    with open("radarNis.csv",'r')as fileReader2:
        radar_Nis = list(csv.reader(fileReader2,delimiter=","))

    radar_Nis_float =  list(map(float,radar_Nis[0]))
    print("num elements:", len(radar_Nis_float))
    radar_95pct = [7.815]*len(radar_Nis_float); 
    ax2 = fig1.add_subplot(122);
    ax2.plot(radar_Nis_float)
    ax2.set_title("Radar NIS");
    ax2.plot(radar_95pct)
    plt.show();


plot_csv();
