import geopandas as gpd  
import geojsonio  
from shapely import geometry  
from matplotlib import pyplot as plt  
import pandas as pd  
import numpy as np  
import glob  
import cv2  
import os  
import json  
import zipfile
# 3 dpf cell plot
egfp_order = glob.glob("./output/egfp_order/*.tif")  
egfp_order.sort()  
egfp_order = egfp_order  
geojsons = glob.glob("./output/sliceGeoJson/*.geojson")  
geojsons.sort()  
geojsons = geojsons  
geojsons
vas_anno = pd.read_csv("../region_annotation_with_feature_F314.csv", index_col = 0)  
vas_anno = vas_anno.drop_duplicates(["region_id"])  
vas_anno = vas_anno[~vas_anno.region_name.isin(['no',np.nan,'?','bottom','Trunk','EYE',"other"])]  
vas_anno.head()
dat = pd.read_csv("../vas_cell_matrix_annotation_F314.csv", index_col = 0)  
dat = dat[~dat["cell_type"].isin(["other", "UnEC"])]  
dat = dat[dat.blood.isin(vas_anno.region_id)]  
dat.shape
dat.replace({"cell_type":{"VEC_1":"VEC","VEC_2":"VEC"}}, inplace=True)
data = dat
data.to_csv("./output/cell_distr/F314_cell_distr.csv", index = False)

for i in range(28):  
    egfp_order_n = egfp_order[i]  
    roi = np.int_(egfp_order_n[-18:-4].split("-")[1:])[::-1]  
    width, height = 3500, 2200  
    start_x, start_y, end_x, end_y = roi[0] - width / 2, roi[1] - height / 2, roi[0] + width / 2, roi[1] + height / 2  
    egfp_image = cv2.imread(egfp_order_n, cv2.IMREAD_GRAYSCALE)  
    plt.imshow(egfp_image * 2.3, cmap="gray")  
    slice = data[(data.x >= start_x) & (data.x <= end_x) & (data.y > start_y) & (data.y < end_y)]  
    colors = ["#ff0000", "#00ff00","#0000ff","#00ffff"]  
    cell_type = ['AEC', 'AngEC', 'CapEC','VEC']  
    for c in range(len(cell_type)):  
        plot_data = slice[slice["cell_type"] == cell_type[c]]  
        plt.scatter(plot_data.x - start_x, plot_data.y - start_y, s=8, label=cell_type[c], color=colors[c])  
        plt.axis("off")  
    plt.tight_layout()  
    plt.axis('off')  
    plt.savefig("./output/cell_distr/{}.tif".format(egfp_order_n[-18:-4]), bbox_inches='tight', pad_inches=0,dpi=300)  
    plt.close()  
    print("done!")

# 6 dpf cell plot
egfp_order = glob.glob("./output/egfp_slices_2/*.tif")  
egfp_order.sort()  
egfp_order
vas_anno = pd.read_csv("../region_annotation_with_feature_F610.csv", index_col = 0)  
vas_anno = vas_anno.drop_duplicates(["region_id"])  
vas_anno = vas_anno[~vas_anno.region_name.isin(['no',np.nan,'?','bottom','Trunk','EYE',"other"])]  
vas_anno.head()
dat = pd.read_csv("../vas_cell_matrix_annotation_F610.csv", index_col = 0)  
dat = dat[~dat["cell_type"].isin(["other", "UnEC"])]  
dat = dat[dat.blood.isin(vas_anno.region_id)]  
dat.shape
dat.replace({"cell_type":{"VEC_1":"VEC","VEC_2":"VEC"}}, inplace=True)
data = dat
data.to_csv("./output/cell_distr/F610_cell_distr.csv", index = False)
for i in range(18):  
    egfp_order_n = egfp_order[i]  
    roi = np.int_(egfp_order_n[-18:-4].split("-")[1:])[::-1]  
    width, height = 2200, 3500  
    start_x, start_y, end_x, end_y = roi[0] - width / 2, roi[1] - height / 2, roi[0] + width / 2, roi[1] + height / 2  
    egfp_image = cv2.imread(egfp_order_n, cv2.IMREAD_GRAYSCALE)  
    plt.imshow(egfp_image * 2.3, cmap="gray")  
    slice = data[(data.x >= start_x) & (data.x <= end_x) & (data.y > start_y) & (data.y < end_y)]  
    colors = ["#ff0000", "#00ff00","#0000ff","#00ffff"]  
    cell_type = ['AEC', 'AngEC', 'CapEC','VEC']  
    for c in range(len(cell_type)):  
        plot_data = slice[slice["cell_type"] == cell_type[c]]  
        plt.scatter(plot_data.x - start_x, plot_data.y - start_y, s=8, label=cell_type[c], color=colors[c])  
        plt.axis("off")  
    plt.tight_layout()  
    plt.axis('off')  
    plt.savefig("./output/cell_distr/{}.tif".format(egfp_order_n[-18:-4]), bbox_inches='tight', pad_inches=0, dpi=300)  
    plt.close()  
    print("done!")

# 11 dpf cell plot
egfp_order = glob.glob("./output/egfp_order/*.tif")  
egfp_order.sort()  
egfp_order
vas_anno = pd.read_csv("../region_annotation_with_feature_F1109.csv", index_col = 0)  
vas_anno = vas_anno.drop_duplicates(["region_id"])  
vas_anno = vas_anno[~vas_anno.region_name.isin(['no',np.nan,'?','bottom','Trunk','EYE',"other"])]  
vas_anno.head()
dat = pd.read_csv("../vas_cell_matrix_annotation_F1109.csv", index_col = 0)  
dat = dat[~dat["cell_type"].isin(["other", "UnEC"])]  
dat = dat[dat.blood.isin(vas_anno.region_id)]  
dat.shape
dat.replace({"cell_type":{"VEC_1":"VEC","VEC_2":"VEC"}}, inplace=True)
data = dat
data.to_csv("./output/cell_distr/F1109_cell_distr.csv", index = False)
for i in range(20):  
    roi = np.int_(egfp_order_n[-18:-4].split("-")[1:])[::-1]  
    width, height = 3500, 2200  
    start_x, start_y, end_x, end_y = roi[0] - width / 2, roi[1] - height / 2, roi[0] + width / 2, roi[1] + height / 2  
    egfp_image = cv2.imread(egfp_order_n, cv2.IMREAD_GRAYSCALE)  
    plt.imshow(egfp_image * 2.3, cmap="gray")  
    slice = data[(data.x >= start_x) & (data.x <= end_x) & (data.y > start_y) & (data.y < end_y)]  
    colors = ["#ff0000", "#00ff00","#0000ff","#00ffff"]  
    cell_type = ['AEC', 'AngEC', 'CapEC','VEC']  
    for c in range(len(cell_type)):  
        plot_data = slice[slice["cell_type"] == cell_type[c]]  
        plt.scatter(plot_data.x - start_x, plot_data.y - start_y, s=8, label=cell_type[c], color=colors[c])  
        plt.axis("off")  
    # ax.legend(['AEC', 'AngEC', 'CapEC', 'VEC'], fontsize=60)  
    plt.tight_layout()  
    plt.axis('off')  
    plt.savefig("./output/cell_distr/{}.tif".format(egfp_order_n[-18:-4]), bbox_inches='tight', pad_inches=0, dpi=300)  
    plt.close()  
    print("done!")
