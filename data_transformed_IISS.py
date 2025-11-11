import json  
import glob  
import numpy as np
import pandas as pd  
import os  
import cv2  
from matplotlib import pyplot as plt  
import random  
import time

# Seperating the vascular counter ROI into a single slice corresponding to the eGFP image.
with open("./emailfiles/11.20250414/Bloods.geojson", "r") as f:  
    x = f.read()  
wgjson = json.loads(x)  
egfp_files = glob.glob("./output/egfp_order/*.tif")  
egfp_files.sort()  
rois = [x[-18:-4] for x in egfp_files]  
for i in range(28):  
    slice_geojson = {'type': 'FeatureCollection','features':[]}  
    filename = "{}.geojson".format(rois[i][:2])  
    roi = np.int_(rois[i].split("-")[1:][::-1])  
    width, height = 3500, 2200  
    start_x, start_y, end_x, end_y = roi[0] - width/2, roi[1] - height/2, roi[0] + width/2, roi[1] + height/2  
    for j in range(len(wgjson["features"])):  
        region = wgjson["features"][j]  
        coord = np.array(region["geometry"]["coordinates"][0])  
        if (all(coord[:,0] > start_x) and all(coord[:,0] < end_x)) and (all(coord[:,1] > start_y) and all(coord[:,1] < end_y)):  
            slice_geojson["features"].append(region)  
    geo = json.dumps(slice_geojson, indent=2)  
    with open("./output/sliceGeoJson/{}".format(filename), "w") as f:  
        f.write(geo)

# Plotting the gene signals on the slices
geneCoord = pd.read_csv("./emailfiles/02.20250313/insitufocus/GeneCoord.csv", index_col = 0)  
egfp_files = glob.glob("./output/egfp_order/*.tif")  
egfp_files.sort()  
rois = [x[-18:-4] for x in egfp_files]
gene_list = geneCoord["feature_name"].unique().tolist()  
gene_list.sort()
color = ["#d32c1f","#acc2d9","#56ae57","#b2996e","#a8ff04","#894585","#d4ffff","#fcfc81","#388004","#efb435","#0c06f7","#3778bf",  
         "#05ffa6","#1f6357","#0cb577","#ff0789","#ff63e9","#430541","#ffb2d0","#ad900d","#6832e3","#850e04","#40fd14","#f6688e",  
         "#76fda8","#014600","#41fdfe","#0c1793","#a50055","#ad03de","#aeff6e","#ff08e8","#fffd01","#0165fc","#f97306","#ffffff",  
         "#0000ff", "#ff0000","#00ff00"]
         
for i in range(28):  
    egfp_file_n = egfp_files[i]  
    roi_str = egfp_file_n[-18:-4]  
    filename = "{}.tif".format(roi_str[:2])  
    roi = np.int_(roi_str.split("-")[1:][::-1])  
    width, height = 3500, 2200  
    start_x, start_y, end_x, end_y = roi[0] - width/2, roi[1] - height/2, roi[0] + width/2, roi[1] + height/2  
    temp = geneCoord[(geneCoord['x_location'] >= start_x) & (geneCoord['x_location'] <= end_x) & (geneCoord['y_location'] >= start_y) & (geneCoord['y_location'] <= end_y)]  
    c = 0  
    pic = cv2.imread(egfp_file_n, cv2.IMREAD_GRAYSCALE)  
    fig = plt.figure(dpi=150, figsize=(width/200, height/200))  
    for gene in gene_list:  
        tmp = temp[temp["feature_name"] == gene]  
        plt.scatter(tmp["x_location"]-start_x, tmp["y_location"]-start_y, s = 8, c = color[c],  
                    alpha = 0.85, linewidths = 0# , label = gene  
                    )  
        c += 1  
    # plt.legend(gene_list)  
    plt.imshow(pic, cmap='gray')  
    plt.axis("off")  
    plt.tight_layout()  
    plt.savefig("./output/allGenePlotOnTheSlices/{}".format(filename), bbox_inches='tight', pad_inches=0)  
    print("{} done! {}".format(filename, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))  
    plt.close(fig)

# The selection of cell located in the vasculature.
vas_cell = pd.read_csv("./emailfiles/11.20250414/cell2blood.csv")  
vas_cell = vas_cell[~vas_cell["blood"].isin(["Non vascular cells", "Perivascular cells"])]  
vas_cell.head()
# plotting the vascular cells
plt.scatter(x = vas_cell["x"], y = vas_cell["y"], color = "red", s=1)  
plt.show()
# checking whether the cells were located on the vascular in a slice
img = cv2.imread("./output/egfp_order/09-14558-06494.tif", cv2.IMREAD_GRAYSCALE)  
vas_cell_t = vas_cell[(vas_cell["x"] > 3000) & (vas_cell["x"] < 10000) & (vas_cell["y"] > 12000) & (vas_cell["y"] < 18000)]  
fig = plt.figure(dpi=150, figsize=(7, 4.4))  
plt.imshow(img, cmap='gray')  
plt.scatter(x = vas_cell_t["x"] - 4750, y = vas_cell_t["y"] - 13450, color = "red", s=1)  
plt.axis('off')  
plt.show()

cell_matrix = pd.read_csv("./emailfiles/04.20250319/cell_matrix.csv", index_col = 0)  
cell_matrix_t = cell_matrix.T
cell_matrix_t_vas_cell = cell_matrix_t[cell_matrix_t.index.isin(vas_cell["cell"])]  
cell_matrix_t_vas_cell.insert(0, "cell", cell_matrix_t_vas_cell.index)
cell_matrix_t_vas_cell = pd.merge(cell_matrix_t_vas_cell, vas_cell, on = "cell")  
cell_matrix_t_vas_cell
