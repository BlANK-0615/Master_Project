import h5py

filename = "/home/s2331261/Master_Project/6_Feature_Selection/params_103.hdf5"

with h5py.File(filename, 'r') as f:
    # 获取并打印文件中的所有顶级组
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]

    # 从文件中获取数据
    data = list(f[a_group_key])

    print(data)
