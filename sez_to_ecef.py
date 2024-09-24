# sez_to_ecef.py
#
# Usage: python3 sez_to_ecef.py o_lat_deg o_lon_deg o_hae_km s_km e_km z_km
#  Converts SEZ to ECEF vector components
# Parameters:
#  o_lat_deg: origin latitude in deg
#  o_lon_deg: origin longitude in deg
#  o_hae_km: origin height above ellipsoid in km
#  s_km: south vector in km
#  e_km: east vector in km
#  z_km: vector normal to the ellipsoid surface in km
#  ...
# Output:
#  A description of the script output XXXXXXXXXXXXXXXXXXXXXXX
#  Inputs an origin location in LLH and a position vector SEZ and prints the ECEF x-component, y-component, and z-component in km
#
# Written by Josefine Krohn
# Other contributors: Brad Denby (format)
#

# import Python modules
import math # math module
import sys # argv

# "constants"
# e.g., R_E_KM = 6378.137

# helper functions
# Matrix Multiplication
def matrix_multiplication(A, B):
    result = [[0 for column_result in range(len(B[0]))] for row_result in range(len(A))] # setting up results vector
    for row_A in range(len(A)): # loop for rows of matrix A
        a = A[row_A] # creating a vector of the current row of A
        for column_B in range(len(B[0])): # loop for columns of matrix B
            b = [row[column_B] for row in B] # creating a vector of the current column of B       
            result_value = 0 # initializing result value
            for column_A in range(len(b)): # loop for # of values in vector b (equals # of values in a)
                value = a[column_A] * b[column_A] # multiplying values in corresponding rows and columns
                result_value = result_value + value # adds the row*column values together
        result[row_A][column_B] = result_value # updates result matrix
    return result  
# Converts LLH to ECEF vector components
def llh_to_ecef(lat_deg, lon_deg, hae_km):
    # import Python modules
    import math # math module
    import sys  # argv
    # "constants"
    R_E_KM = 6378.1363
    E_E = 0.081819221456
    # converting angles in deg to rad
    lat_rad = lat_deg*math.pi/180
    lon_rad = lon_deg*math.pi/180
    # calculations
    denom = math.sqrt(1 - (E_E**2)*(math.sin(lat_rad)**2))
    c_E = R_E_KM/denom
    s_E = R_E_KM*(1 - E_E**2)/denom
    r_x_km = (c_E + hae_km)*math.cos(lat_rad)*math.cos(lon_rad)
    r_y_km = (c_E + hae_km)*math.cos(lat_rad)*math.sin(lon_rad)
    r_z_km = (s_E + hae_km)*math.sin(lat_rad)
    r_vec = [[r_x_km], [r_y_km], [r_z_km]]
    return r_vec # ECEF x-component, y-component, and z-component in km

# initialize script arguments
o_lat_deg = float('nan') # origin latitude in deg
o_lon_deg = float('nan') # origin longitude in deg
o_hae_km = float('nan') # origin height above ellipsoid in km
s_km = float('nan') # south vector in km
e_km = float('nan') # east vector in km
z_km = float('nan') # vector normal to the ellipsoid surface in km

# parse script arguments
if len(sys.argv)==7:
  o_lat_deg = float(sys.argv[1])
  o_lon_deg = float(sys.argv[2])
  o_hae_km = float(sys.argv[3])
  s_km = float(sys.argv[4])
  e_km = float(sys.argv[5])
  z_km = float(sys.argv[6])
else:
  print(\
   'Usage: '\
   'python3 sez_to_ecef.py o_lat_deg o_lon_deg o_hae_km s_km e_km z_km'\
  )
  exit()


# write script below this line

o_lat_rad = o_lat_deg*math.pi/180 # converting from degrees to radians
o_lon_rad = o_lon_deg*math.pi/180 # converting from degrees to radians

r_sez = [[s_km], [e_km], [z_km]] # SEZ vector


# SEZ to ECEF rotation
phi = o_lat_rad
theta = o_lon_rad
R_z = [[math.cos(theta), -math.sin(theta), 0], [math.sin(theta), math.cos(theta), 0], [0, 0, 1]]
R_y = [[math.sin(phi), 0, math.cos(phi)], [0, 1, 0], [-math.cos(phi), 0, math.sin(phi)]]
R_y_r_sez = matrix_multiplication(R_y,r_sez) # R_y(90deg - phi)*r_sez
r_ECEF = matrix_multiplication(R_z,R_y_r_sez) # r_ECEF = R_z(theta)*R_y(90deg - phi)*r_sez


o_r_ECEF = llh_to_ecef(o_lat_deg, o_lon_deg, o_hae_km) # converting origin from LLH to ECEF

#add ecef vector to sez origin
ecef_x_km = o_r_ECEF[0][0]+r_ECEF[0][0]
ecef_y_km = o_r_ECEF[1][0]+r_ECEF[1][0]
ecef_z_km = o_r_ECEF[2][0]+r_ECEF[2][0]

print(ecef_x_km)
print(ecef_y_km)
print(ecef_z_km)