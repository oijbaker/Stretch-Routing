from pandas import read_csv
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def periodic_dist(p1, p2, box_dims):
    dx = np.abs(p1[0] - p2[0])
    dx = np.minimum(dx, box_dims[0] - dx)
    dy = np.abs(p1[1] - p2[1])
    dy = np.minimum(dy, box_dims[1] - dy)
    return np.sqrt(dx**2 + dy**2)

def draw_toroidal_edge(ax, p1, p2, box_dims, **kwargs):
    x1, y1 = p1
    x2, y2 = p2
    a, b = box_dims
    
    dx = x2 - x1
    dy = y2 - y1
    wrap_x = abs(dx) > a / 2
    wrap_y = abs(dy) > b / 2
    
    x2_col = x2 - a if (wrap_x and dx > 0) else (x2 + a if (wrap_x and dx < 0) else x2)
    x1_col = x1 + a if (wrap_x and dx > 0) else (x1 - a if (wrap_x and dx < 0) else x1)
    y2_col = y2 - b if (wrap_y and dy > 0) else (y2 + b if (wrap_y and dy < 0) else y2)
    y1_col = y1 + b if (wrap_y and dy > 0) else (y1 - b if (wrap_y and dy < 0) else y1)
    
    if not wrap_x and not wrap_y:
        ax.plot([x1, x2], [y1, y2], **kwargs)
    else:
        ax.plot([x1, x2_col], [y1, y2_col], **kwargs)
        ax.plot([x1_col, x2], [y1_col, y2], **kwargs)

def get_hex_neighbors(r, c, M, N):
    if r % 2 == 0:
        neighbors = [(r, c-1), (r, c+1), (r-1, c-1), (r-1, c), (r+1, c-1), (r+1, c)]
    else:
        neighbors = [(r, c-1), (r, c+1), (r-1, c), (r-1, c+1), (r+1, c), (r+1, c+1)]
    return [((nr % N), (nc % M)) for nr, nc in neighbors]

def get_hex_vertices(center, L):
    vertices = []
    for i in range(6):
        angle = np.pi / 6 + i * np.pi / 3
        vertices.append((center[0] + L * np.cos(angle), center[1] + L * np.sin(angle)))
    return vertices


def plot_with_hexagons(G, G_sub, points, a, b):
    box_dims = np.array([a,b])
    fig, ax = plt.subplots(1, 2, figsize=(20, 10 * (b/a)))
    ax[0].set_xlim(0, a)
    ax[0].set_ylim(0, b)
    
    for edge in G_sub.edges():
        draw_toroidal_edge(ax[0], points[edge[0]], points[edge[1]], box_dims, 
                        color='royalblue', alpha=0.7, lw=1)


    # Draw the hexagonal grid
    for r in range(N):
        for c in range(M):
            y = r * 1.5 * L
            x = c * np.sqrt(3.0) * L + (0.5 * np.sqrt(3.0) * L if r % 2 == 1 else 0.0)
            center = (x, y)
            vertices = get_hex_vertices(center, L)
            hexagon = plt.Polygon(vertices, closed=True, edgecolor='gray', facecolor='none', alpha=0.3)
            ax[0].add_patch(hexagon)
        
    ax[0].scatter(points[:,0], points[:,1], s=10)

    ax[0].set_title(f"Compressed Graph")
    ax[0].legend(loc='upper right')
    ax[0].set_aspect('equal')

    # Draw the base graph in the second subplot
    for node in G.nodes():
        for neighbor in G.neighbors(node):
            if node < neighbor:  # Avoid double plotting
                draw_toroidal_edge(ax[1], points[node], points[neighbor], box_dims, 
                                color='lightgray', alpha=0.7, lw=1)
                
    plt.scatter(points[:, 0], points[:, 1], c='black', s=5, zorder=2)
                
    ax[1].set_title("Base Graph (SRGG)", fontsize=14)
    ax[1].set_xlim(0, a)
    ax[1].set_ylim(0, b)
    ax[1].set_aspect('equal')

    plt.grid(True, linestyle=':', alpha=0.3)
    plt.show()




def plot(G, G_sub, points, a, b):
    box_dims = np.array([a,b])
    fig, ax = plt.subplots(1, 2, figsize=(20, 10 * (b/a)))
    ax[0].set_xlim(0, a)
    ax[0].set_ylim(0, b)
    
    for edge in G_sub.edges():
        draw_toroidal_edge(ax[0], points[edge[0]], points[edge[1]], box_dims, 
                        color='royalblue', alpha=0.7, lw=1)
        
    ax[0].scatter(points[:,0], points[:,1], s=10)

    ax[0].set_title(f"Compressed Graph")
    ax[0].legend(loc='upper right')
    ax[0].set_aspect('equal')

    # Draw the base graph in the second subplot
    for node in G.nodes():
        for neighbor in G.neighbors(node):
            if node < neighbor:  # Avoid double plotting
                draw_toroidal_edge(ax[1], points[node], points[neighbor], box_dims, 
                                color='lightgray', alpha=0.7, lw=1)
                
    plt.scatter(points[:, 0], points[:, 1], c='black', s=5, zorder=2)
                
    ax[1].set_title("Base Graph (SRGG)", fontsize=14)
    ax[1].set_xlim(0, a)
    ax[1].set_ylim(0, b)
    ax[1].set_aspect('equal')

    plt.grid(True, linestyle=':', alpha=0.3)
    plt.show()


print("reading")
dfG = read_csv("output/graph_edges.csv", names=["0", "1"])
dfG_sub = read_csv("output/subgraph_edges.csv", names=["0", "1"])
points = read_csv("output/points.csv", names=["x", "y"])

G = np.array(dfG)
G_sub = np.array(dfG_sub)
points = np.array(points)

G = nx.from_edgelist(G)
G_sub = nx.from_edgelist(G_sub)

L = 0.11
a = 1.0  # Fix width to 1

# Solve for M (number of columns) and snap to closest integer
M = int(np.round(1.0 / (L * np.sqrt(3))))
M = max(1, M) # Ensure at least 1 column

# Re-adjust L slightly so M hexagons fit perfectly across width 1.0
L = 1.0 / (M * np.sqrt(3))

# Solve for N (number of rows) to get height b close to 1.0
# b = N * 1.5 * L -> N = 1.0 / (1.5 * L)
ideal_N = 1.0 / (1.5 * L)
# N must be an even integer for vertical periodic stitching
N = int(2 * np.round(ideal_N / 2))
N = max(1, N)

# Final height of the torus
b = N * 1.5 * L

plot_with_hexagons(G, G_sub, points, a, b)

