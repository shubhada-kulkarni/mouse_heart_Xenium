# python script to remove bad cell from the segmentation that was causing
# yellowing of whole image in Xenium Explorer
#
import sys, os
import geojson

cellfile = sys.argv[1] #"cell_boundaries_18S.geojson"
nucleusfile = sys.argv[2] #"nuclei_DAPI.geojson"

outcellfile = os.path.dirname(cellfile) + "/cell_boundaries_18S_mod.geojson"
outnucleusfile = os.path.dirname(cellfile) + "/nuclei_DAPI_mod.geojson"

#################################################################
## cell
# file reading
with open(cellfile) as fcell:
    gcell = geojson.load(fcell)

# pop two times - first for small round circle, second for entire slide polygon
gcell['features'].pop(0)
gcell['features'].pop(0)

# write modified object
with open(outcellfile, 'w') as file:
    geojson.dump(gcell, file)

#################################################################
## nucleus
# file reading
with open(nucleusfile) as fnuc:
    gnuc = geojson.load(fnuc)

# pop two times - first for small round circle, second for entire slide polygon
gnuc['features'].pop(0)
gnuc['features'].pop(0)

# write modified object
with open(outnucleusfile, 'w') as file:
    geojson.dump(gnuc, file)



'''
fin = "/prj/XeniumProbeDesign/mouse_heart_03012025/output-XETG00046__0046388__Region_1__20241219__110827/cell_boundaries_18S.geojson" # for cells
# fin = "/prj/XeniumProbeDesign/mouse_heart_03012025/output-XETG00046__0046388__Region_1__20241219__110827/nuclei_DAPI.geojson"   # for nuclei

with open(fin) as f:
    gj = geojson.load(f)

# for nuclei
for i in range(2, len(gj['features'])):
    cellid = gj['features'][i]["id"]
    area = gj['features'][i]['properties']['measurements']["Nucleus: Area"]
    mean_18s = gj['features'][i]['properties']['measurements']["Nucleus: 18S mean"]
    max_18s = gj['features'][i]['properties']['measurements']['Nucleus: 18S max']
    print(i, cellid, area, mean_18s, max_18s)


# for cells
for i in range(2, len(gj['features'])):
    cellid = gj['features'][i]["id"]
    area = gj['features'][i]['properties']['measurements']["Cell: Area"]
    mean_18s = gj['features'][i]['properties']['measurements']["Cell: 18S mean"]
    max_18s = gj['features'][i]['properties']['measurements']['Cell: 18S max']
    
    x = []  # for coord info storage
    y = []  
    coord = gj['features'][i]["geometry"]["coordinates"][0]
    for each in coord:
        x.append(each[0])
        y.append(each[1])

    print(i, cellid, area, mean_18s, max_18s, min(x), max(x), min(y), max(y))
'''
