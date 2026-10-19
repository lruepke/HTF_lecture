import os

def write_stl(vti_file, stl_file):
    # 1. load vti file
    data            = OpenDataFile(vti_file)
    # 2. clip at some intermediate value (we have 0 and 255 as pores and grains)
    clip1           = Clip(data, ClipType = 'Scalar', Scalars = ['CELLS', 'im'], Value = 125, Invert = 1)
    # 3. make a surface of the remaining grains
    extractSurface1 = ExtractSurface(clip1)
    # 4. and triangulate it for stl export
    triangulate1    = Triangulate(extractSurface1)

    # 5. finally save it as an stl file
    SaveData(stl_file, proxy = triangulate1)

# main part


print(f"Current working directory: {os.getcwd()}")
vti_file = 'porous_model_spheres.vti'      # input .vti file
stl_file = '../constant/triSurface/porous_model_spheres.stl'  # output .stl file
# call function
write_stl(vti_file, stl_file)