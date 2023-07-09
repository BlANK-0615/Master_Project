import csv

with open('/home/s2331261/Master_Project/6_Feature_Selection/CC.txt', 'r') as in_file:
    stripped = (line.strip() for line in in_file)
    lines = (line.split(":") for line in stripped if line)
    with open('CC.csv', 'w') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(('file name', 'conformations count'))  # 写入列名
        writer.writerows(lines) 




