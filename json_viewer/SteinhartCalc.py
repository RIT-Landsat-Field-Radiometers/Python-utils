# import numpy as np
#
#
# def main():
#     points = [
#         (0, 94980.00),
#         (25, 30000.00),
#         (50, 10969.00),
#         (70, 5357.40)]
#
#     mat1 = np.array([[np.log(points[x][1]) ** y for y in range(4)] for x in range(len(points))])
#     inv = np.linalg.inv(mat1)
#     mat2 = np.array([[1.0 / (273.15 + points[x][0])] for x in range(len(points))])
#     coeff = np.matmul(inv, mat2)
#
#     print("Coefficients")
#     names = ['A', 'B', 'C', 'D']
#     for i in range(len(names)):
#         print(names[i], end="=")
#         print(coeff[i][0])
#     print("Excel (Kelvin)")
#     print("=1/({A} + {B} * LN(A1) + {C} * POW(LN(A1),2) + {D} * POW(LN(A1),3))".format(A=coeff[0][0], B=coeff[1][0], C=coeff[2][0], D=coeff[3][0]))
#     # print(coeff)


if __name__ == "__main__":
    with open('./log.txt', 'r') as ifp:
        with open('./out.csv', 'w') as ofp:
            line = ifp.readline()
            while line:
                parts = line.strip().split(' ')
                parts = [parts[0].split('.')[0].split(':')[-1], parts[2], parts[-1]]
                print(",".join(parts))
                ofp.write(",".join(parts) + "\n")
                line = ifp.readline()


    # with open("/media/reedt/XTRM-Q/Data/work/DIRS/Temp/comp/1636063972011.gz", 'rb') as ifp:
    #     odata = ifp.read()
    #     with open("/media/reedt/XTRM-Q/Data/work/DIRS/Temp/comp/prepend.gz", 'wb') as ofp:
    #         ofp.write(bytearray([0 for _ in range(512)]))
    #         ofp.write(odata)

    # main()
