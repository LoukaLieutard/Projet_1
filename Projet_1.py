import random
import math



def xy_to_ij(lst_xy, n_size_matrix):
    # Transforme les coordonnée réelles en coordonnées de matrice
    j = int(lst_xy[0])
    i = n_size_matrix - 1 - int(lst_xy[1])
    return[i, j]
print(xy_to_ij([2.5, 2.5], 5))

def calculate_fraction_E_in_cells(lst_coords_in_cell, matE, Xwidth, n_cells) :
    mat_size = 5
    if len(matE) != mat_size :
        print("\n\n*** ATTENTION ***\nLa taille de la matrice fournie en argument de calculate_fraction_E_in_cells n'a pas l'air correcte.\nSortie anormale de la fonction.\n\n\n")
        return None
    if Xwidth != 100 :
        print("\n\n*** ATTENTION ***\nLa largeur de detecteur recue est",Xwidth,"cm au lieu de 100 cm.\nSortie anormale de la fonction.\n\n\n")
        return None
    if not (5 <= n_cells <= 50) :
        print("\n\n*** ATTENTION ***\nLe nombre", n_cells," de cellules recu est invalide.\nSortie anormale de la fonction.\n\n\n")
        return None
    #Quick trick to see an actual center for large n_cells and still keep fluctuations for low n_cells:
    #nbTirages = int(100 * (n_cells / 20))
    #Let's not use it eventually, fluctuations may be closer to reality this way.
    nbTirages = 100
    #Calo config rather than pixel.
    #In this situation, shower width = 2.5 cm (a-la Moliere radius).
    #In variable sigma though, shower width is stored in units of cell width.
    #Therefore sigma = 2.5cm / cell_width, ie sigma = 2.5cm * n_cells / detector_width.
    #1st limit is when 3 sigmas (in cm) becomes larger than 5x5 matrix diagonal = 2.5 cell_width sqrt(2). This gives n_cells < 2.5 Xwidth sqrt(2) / 3*2.5cm ~ 47.
    #2nd limit is when 3 sigmas (in cm) is all in a single tower, ie n_cells > Xwidth / 3sigma ~ 13.
    sigma = 2.5 * n_cells / Xwidth
    #Just in case the user doesn't provide an empty matrix...
    matRdm = create_empty_matrix(mat_size)
    for i in range (0,nbTirages) :
        xRdm = random.gauss(float(mat_size//2) + lst_coords_in_cell[0], sigma)
        yRdm = random.gauss(float(mat_size//2) + lst_coords_in_cell[1], sigma)
        #In rare cases the random number is picked outside the interval [0,mat_size[, which triggers an out-of-range error when used in xy_to_ij / matRdm.
        if xRdm < 0 :
           xRdm = 0
        if yRdm < 0 :
           yRdm = 0
        if xRdm >= mat_size :
           xRdm = mat_size-0.001
        if yRdm >= mat_size :
           yRdm = mat_size-0.001
        lst_Rdm = xy_to_ij([int(xRdm),int(yRdm)],mat_size)
        matRdm[lst_Rdm[0]][lst_Rdm[1]] += 1
    for i in range(mat_size) :
        for j in range(mat_size) :
            matE[i][j] = matRdm[i][j] / nbTirages
    return None

def create_empty_matrix(nbCells, marge = 0):
    # Crée une matrice vide avec une taille et une taille de marge définie
    matrix = []
    for i in range(nbCells+2*marge):
        matrix.append([])
        for j in range(nbCells+2*marge):
            matrix[i].append(0)
    return matrix

def print_matrix(matrix):
    # Affiche le contenu de la matrice de taille quelconque
    for i in matrix:
        for j in i:
            print(int(j*100)/100, end = " ")
        print()

def print_matrix_bool(matrix):
    # Affiche un point si une cellule = 0 et un # si non
    for i in matrix:
        for j in i:
            if j == 0:
                print(".", end = " ")
            else:
                print("#", end = " ")
        print()

def print_matrix_energy(matrix):
    # Affiche différents caracteres en fonction de l'energie dans la cellule
    for i in matrix:
        for j in i:
            if j > 0.75:
                print("#", end = "  ")
            elif j > 0.5:
                print("*", end = "  ")
            elif j > 0.25:
                print("o", end = "  ")
            elif j > 0.01:
                print(".", end = "  ")
            else:
                print("--", end = " ")       
        print()

def generate_partic_coords(Xwidth):
    # Génere une particule et renvoie ses coordonnées
    return [random.random()*Xwidth, random.random()*Xwidth]

def calc_all_coords(lst_coords_in_det, cell_size):
    # Renvoie les coordonées de la cellule touchée et les coordonée a l'interieur de la cellule
    coord_det = lst_coords_in_det.copy()
    cell_coord = [0,0]
    coord_in_cell = [0,0]
    while coord_det[0] > 0:
        coord_det[0] -= cell_size
        cell_coord[0] += 1
    coord_in_cell[0] = (coord_det[0]+cell_size)/cell_size
    while coord_det[1] > 0:
        coord_det[1] -= cell_size
        cell_coord[1] += 1
    coord_in_cell[1] = coord_det[1]+cell_size
    return [cell_coord, coord_in_cell]

def smear_energy(energy_matrix):
    matrix_copy = energy_matrix.copy()
    for i in matrix_copy:
        for j in i:
            j *= random.gauss(1, 0.05)
    return matrix_copy

def apply_threshold(energy_matrix):
    for i in energy_matrix:
        for j in i:
            if j < 0.05:
                j = 0.0

def simulate_signal(lst_coords_in_cell, Xwidth, n_cells):
    #crée une matrice 5x5 vide puis simule une particule qui vient la percuter en son centre et aplique les fonction smear_energy et apply_threshold a la matrice resultant 
    matE =  []
    for i in range(5):
        matE.append([])
        for j in range(5):
            matE[i].append(0)
    print(lst_coords_in_cell)
    print(Xwidth)
    print(n_cells)
    calculate_fraction_E_in_cells(lst_coords_in_cell, matE, Xwidth, n_cells)
    print(matE)
    matE = smear_energy(matE)
    apply_threshold(matE)
    for i in matE :
        print(i)
        print()
    return matE

def launch_particle_on_detector(Xwidth, N, n, matrix_N_2n):
    cell_size = Xwidth / N
    lst_temp = generate_partic_coords(Xwidth)
    x_y = xy_to_ij(calc_all_coords(lst_temp, cell_size)[0])
    a_b = calc_all_coords(lst_temp, cell_size)[1]
    a_b[1] = 1 - a_b[1]
    matrix_E = simulate_signal(a_b, Xwidth, N**2)

#simulate_signal([0.5, 0.5], 100, 20)