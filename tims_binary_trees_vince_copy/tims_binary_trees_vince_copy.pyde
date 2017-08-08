def setup():
    size(600, 600, P2D)
    
    global n_ary; # How many branches at each node.
    n_ary = 2;
    
    global canopies; # Stores the nodes as the various levels.
    canopies = [[] for i in range(20)] # max depth of 20
    
    # These keep track of the size of the image in user space.
    global min_x, max_x, min_y, max_y;
    import math
    
    
# Parameters for tree:
# 1 basis             :  The list of affine transformations (the basis).
# 2 trunk             :  The trunk of the tree.
# 3 depth             :  The number of levels to create.
# 4 bgcolor           :  Background color.
# 5 init_circle_radius:  How big the first nodes are.
# 6 radius_ratio      :  How the sizes of the nodes change with each level.
# 7 node_colors       :  List of the colors of the nodes at each level.    
# 8 init_branch_width :  How wide the trunk is.
# 9 branch_ratio      :  How the width of the branches change with each level.
#10 branch_colors     :  List of the colors of the branches at each level.

def draw():
    if frameCount < 2:
        tree(
              [R(25), IN([1,2], [1.5,-1], 0.5,0.75)], # 1 basis
              [0,1], # 2 trunk
              6,  # 3 depth
              color(50,50,50),  # 4 bgcolor
              3,  # 5 init_circle_radius
              0.9, # 6 radius_ratio
              # 7 node_colors
              [color(255,100,0,255 - 15 * alpha) for alpha in range(20)], 
              3, # 8 init_branch_width
              0.85, # 9 branch_ratio
              # 10 branch_colors
              #[color(255,255,0, 255 - 15 * alpha) for alpha in range(20)] 
              [color(0,255,255), color(255,0,255), color(255,255,0), 
               color(0,0,255), color(0,255,0), color(255,0,0)]
              );
        saveFrame("frames//####.png")

    else:
        noLoop();

# Affine transformations are represented as a list [a, b, c, d, e, f].
def apply (aff, vec):
    #print aff, vec
    return [aff[0] * vec[0] + aff[1] * vec[1] + aff[4], 
        aff[2] * vec[0] + aff[3] * vec[1] + aff[5]]
    
# Affine function composition.    
def compose (aff_1, aff_2):
    return [aff_1[0] * aff_2[0] + aff_1[1] * aff_2[2], 
              aff_1[0] * aff_2[1] + aff_1[1] * aff_2[3], 
              aff_1[2] * aff_2[0] + aff_1[3] * aff_2[2], 
              aff_1[2] * aff_2[1] + aff_1[3] * aff_2[3], 
              aff_1[0] * aff_2[4] + aff_1[1] * aff_2[5] + aff_1[4],
              aff_1[2] * aff_2[4] + aff_1[3] * aff_2[5] + aff_1[5]]
    
# Now defining the various ways to specify transformations.  

# Defines a rotation matrix given an angle in degrees.
def R(deg):      
    theta = radians(deg);
    return [cos(theta), -sin(theta), sin(theta), cos(theta), 0, 0];

def transpose(aff):
    return [aff[0], aff[2], aff[1], aff[3], aff[4], aff[5]];

def inverse(aff):
    det = float(aff[0] * aff[3] - aff[1] * aff[2]); # Use "float" in the event aff has integer components.
    linear = [aff[3]/det, -aff[1]/det, -aff[2]/det, aff[0]/det];
    return linear + apply(linear + [0,0], [aff[4], aff[5]]);

# Now defining the four types of matrices (excludes the zero matrix). 
# Singular and nondefective.
def SN(evec_1, evec_2, eval_1):
    P = transpose(evec_1 + evec_2 + [0, 0]);
    return compose(compose(P, [eval_1, 0, 0, 0, 0, 0]), inverse(P));

# Singular and defective.
def SD(evec_1, evec_2):
    P = transpose(evec_1 + evec_2 + [0, 0]);
    return compose(compose(P, [0, 1, 0, 0, 0, 0]), inverse(P));

# Invertible and nondefective.  Recall that eigenvalues/vectors may be complex conjugates.
def get_real(aff):
    return [aff[0].real, aff[1].real, aff[2].real, aff[3].real, aff[4].real, aff[5].real];

def IN(evec_1, evec_2, eval_1, eval_2):
    P = transpose(evec_1 + evec_2 + [0, 0]);
    return compose(compose(P, [eval_1, 0, 0, eval_2, 0, 0]), inverse(P));
    
# Computes the next node in a path.     
def next (node, aff, trunk):
    global min_x, max_x, min_y, max_y;

    # Computes the next branch.
    # A node is a list with this structure: [[0, 1], [2, 3, 4, 5, 6, 7]].
    # The first element is the coordinates of the node, the second is the direction (eta).
    eta = compose(aff, node[1]);
    branch = apply(eta, trunk[0])
    # Performs vector addition.
    nextnode = [node[0][0] + branch[0], node[0][1] + branch[1]];
    
    # updates range of x and y values.
    min_x = min(min_x, nextnode[0]);
    max_x = max(max_x, nextnode[0]);
    min_y = min(min_y, nextnode[1]);
    max_y = max(max_y, nextnode[1]);
    
    return [nextnode, eta] # The next node.
    
def myline (x_1, y_1, x_2, y_2):
    global  min_x, max_x, min_y, max_y;
    padding = 0.05; # How much the image is padded as a fraction of the width/height.
    
    # These quantities are used frequently.
    dx = max_x - min_x;
    dy = max_y - min_y;
    aspect_ratio = dx / dy;

    if aspect_ratio > float(width)/height:
        # Larger aspect ratio than the screen; need y padding.
        # First calculate the width and height.
        x_width = (1 - 2 * padding) * width;
        y_height = x_width / aspect_ratio;
        
        # Now find the borders.
        y_border = (height - y_height) / 2.
        x_border = padding * width;
                
        line(x_border + (x_1 - min_x) / dx * x_width,
            height - (y_border + (y_1 - min_y) / dy * y_height),
            x_border + (x_2 - min_x) / dx * x_width,
            height - (y_border + (y_2 - min_y) / dy * y_height));
    else:
        # Smaller aspect ratio than the screen; need x padding.
        # First calculate the height and width.
        y_height = (1 - 2 * padding) * height;
        x_width =  aspect_ratio * y_height;
       
        # Now find the borders.
        x_border = (width - x_width) / 2;
        y_border = padding * height;
        
        line(x_border + (x_1 - min_x) / dx * x_width,
            height - (y_border + (y_1 - min_y) / dy * y_height),
            x_border + (x_2 - min_x) / dx * x_width,
            height - (y_border + (y_2 - min_y) / dy * y_height));
            
def mycircle (x_1, y_1, r):
    global  min_x, max_x, min_y, max_y;
    padding = 0.05; # How much the image is padded as a fraction of the width/height.
    
    # These quantities are used frequently.
    dx = max_x - min_x;
    dy = max_y - min_y;
    aspect_ratio = dx / dy;

    if aspect_ratio > float(width)/height:
        # Larger aspect ratio than the screen; need y padding
        # First calculate the width and height.
        x_width = (1 - 2 * padding) * width;
        y_height =  x_width / aspect_ratio;
        
        # Now find the borders.
        y_border = (height - y_height) / 2.
        x_border = padding * width;
        
        ellipse(x_border + (x_1 - min_x) / dx * x_width,
            height - (y_border + (y_1 - min_y) / dy * y_height),
            r, r);
    else:
        # Smaller aspect ratio than the screen; need x padding
        # First calculate the height and width.
        y_height = (1 -2 * padding) * height;
        x_width =  aspect_ratio * y_height;
       
        # Now find the borders.
        x_border = (width - x_width) / 2;
        y_border = padding * height;
        
        ellipse(x_border + (x_1 - min_x) / dx * x_width,
            height - (y_border + (y_1 - min_y) / dy * y_height),
            r, r);
# Parameters for tree:
# 1 basis             :  The list of affine transformations (the basis).
# 2 trunk             :  The trunk of the tree.
# 3 depth             :  The number of levels to create.
# 4 bgcolor           :  Background color.
# 5 init_circle_radius:  How big the first nodes are.
# 6 radius_ratio      :  How the sizes of the nodes change with each level.
# 7 node_colors       :  List of the colors of the nodes at each level.    
# 8 init_branch_width :  How wide the trunk is.
# 9 branch_ratio      :  How the width of the branches change with each level.
#10 branch_colors     :  List of the colors of the branches at each level.
                        
def tree (basis, trunk, depth, bgcolor,
          init_circle_radius, radius_ratio, node_colors,
          init_branch_width, branch_ratio, branch_colors): 
    
    global n_ary, canopies, min_x, max_x, min_y, max_y;
    
    background(bgcolor);
    
    # Initializes the bounding box parameters.
    min_x = min(0,trunk[0]);
    max_x = max(0,trunk[0]);
    min_y = min(0,trunk[1]);
    max_y = max(0,trunk[1]);
    
    # Creates the trunk.
    trunk = [trunk, [1,0,0,1,0,0]];
    canopies[0] = [trunk];
    
    # Now creates the nodes at various depths.
    for level in range(depth):
        canopies[level + 1] = [next(node, basis[i], trunk) for node in canopies[level]  for i in range(n_ary)];
    
    # When doing linear interpolation, the tree might lie on a line.
    # We must avoid dividing by 0 in the circle/line procedures.
    # This happens rarely, so it's not worth trying to avoid this cleverly.
    if min_x == max_x:
        max_x = min_x + 0.0000001;
    if min_y == max_y:
        max_y = min_y + 0.0000001;
    
    # Draws the trunk.
    strokeCap(ROUND);
    stroke(color(branch_colors[0]));
    strokeWeight(init_branch_width);
    myline(canopies[0][0][0][0], canopies[0][0][0][1],0,0);
    
    # Drawing the branches.            
    for level in range(depth):
        # Sets branch color and width.
        stroke(color(branch_colors[level]));
        strokeWeight(init_branch_width * (branch_ratio ** level));
        # Draws the branches.
        [myline(canopies[level + 1][n][0][0], canopies[level + 1][n][0][1], 
                canopies[level][n/n_ary][0][0], canopies[level][n/n_ary][0][1]) for n in range(n_ary**(level+1))]
        
    # Drawing the nodes.   
    for level in range(depth):
        # Sets node color and size, then renders them.
        strokeWeight(0);
        stroke(color(node_colors[level]));
        fill(color(node_colors[level]));
        [mycircle(canopies[level + 1][n][0][0], canopies[level + 1][n][0][1],
                  init_circle_radius * (radius_ratio ** level) ) for n in range(n_ary ** (level+1))];


    