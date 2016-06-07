This folder contains data used for experimets in the following paper:

Moral Lineage Tracing. F. Jug*, E. Levinkov*, C. Blasse, E. Myers, and B. Andres
In: IEEE International Conference on Computer Vision and Pattern Recognition, 2016

*contributed equally

Each dataset consists of two files describing nodes and edges of the graph.


nodes.csv constains information about each vertex per line in the following format:

frame_number node_id_in_frame x y p
frame_number node_id_in_frame x y p
frame_number node_id_in_frame x y p
frame_number node_id_in_frame x y p

where x and y are coordinates of the center of mass of the supervoxel in frame and p is the probability of the superpixel to be born or terminated in the current frame. The latter is estimated siply based on the distance from the image border. In fact, in our experiments reported in the paper we used default weights for any node, so this field is not really used in the current system.


edges.csv contains information about each edge per line in the following format:

frame_number_0 node_id_0 frame_number_1 node_id_1 p
frame_number_0 node_id_0 frame_number_1 node_id_1 p
frame_number_0 node_id_0 frame_number_1 node_id_1 p
frame_number_0 node_id_0 frame_number_1 node_id_1 p

frame_number and node_id strictly coincides with te indexing defined in nodes.csv. Here, p is the estimated probability of an edge being cut.

To easily load the data, just use the code in lineage/problem.hxx