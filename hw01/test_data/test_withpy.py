# 定义文件路径
#file_path = "C:/Users/Acer/Documents/GitHub/Geo1004/Geo1004/hw01/NL.IMBAG.Pand.0503100000020110-0.obj"
file_path = "/Users/Acer/Documents/GitHub/Geo1004/Geo1004/hw01/NL.IMBAG.Pand.0503100000020110-0.obj"
# 使用 with 语句打开文件，确保文件正确关闭
with open(file_path, 'r') as file:
    # 逐行读取文件内容
    for line in file:
        print(line.strip())  # 打印每一行，并去除两端的空白字符
