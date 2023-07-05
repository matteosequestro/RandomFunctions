clear all; clc; close all;

path        = 'C:\Users\Social-Lab\Desktop\Matteo\AAA_Studies\2022_July-Oct_Elevator\data\AAC_Exp2\data\behavioral';
folders     = dir([path, '\Elevator*']);

outpath     = 'C:\Users\Social-Lab\Desktop\Matteo\AAA_Studies\2022_July-Oct_Elevator\data\AAC_Exp2\data';
mkdir([outpath, '\headCorrected2'])

load("C:\Users\Social-Lab\Desktop\Matteo\AAA_Studies\2022_July-Oct_Elevator\properpaths_headtracking.mat");
properPaths = this; clear this;

for version = 1 height(folders)
    
    version =1;
    IdsPath = [path, '\' folders(version).name];
    verIds = dir([IdsPath, '\VR*']);
    
    for subject = 1:height(verIds)
        subject = 1;
        trackers = dir([IdsPath ,'\' verIds(subject).name, '\S001\trackers']);




    end









end



