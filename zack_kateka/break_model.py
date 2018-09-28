from change_biotapestry import remove_biotapestry


"""
This file will remove connections from your network. Here, remove_connections is a list of tuples.
Each tuple is of the form (Source, Target), and represents a connection to remove.
For example, if remove_connections = [(1,2), (5,6)], the connection starting at 1 and ending at 2,
and the connection starting at 5 and ending at 6 will both be removed.
This will not alter your working biotapestry CSV, but will generate a new CSV for the broken model in your working directory.
Run this script when you set remove_connections.
"""

#CHANGE ME
remove_connections = [(4,1), (4,6), (4,2), (6,4),(3,4),(7,7)]

#do NOT change this
model_name = "pathway"

remove_biotapestry(remove_connections, model_name + "_biotapestry.csv", model_name + "_broken.csv")
