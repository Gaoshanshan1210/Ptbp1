import pandas as pd
import argparse
import random
import os

class iTOL_annotation:
    def __init__(self, fi, fo, sep):
        self.fi = fi
        self.fo = fo
        self.sep = sep
        self.annotation_prefix = """DATASET_COLORSTRIP
#In colored strip datasets, each ID is associated to a color box/strip and can have an optional label. Color can be specified in hexadecimal, RGB or RGBA notation. When using RGB or RGBA notation, you cannot use COMMA as the dataset separator

#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.

#SEPARATOR TAB
#SEPARATOR COMMA
SEPARATOR SPACE

#label is used in the legend table (can be changed later)
DATASET_LABEL color_strip1

#dataset color (can be changed later)
COLOR #ff0000

#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#

#If COLOR_BRANCHES is set to 1, branches of the tree will be colored according to the colors of the strips above the leaves.
#When all children of a node have the same color, it will be colored the same, ie. the color will propagate inwards towards the root.
COLOR_BRANCHES 0


#=================================================================#
#     all other optional settings can be set or changed later     #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#

#Each dataset can have a legend, which is defined using LEGEND_XXX fields below
#For each row in the legend, there should be one shape, color and label.
#Optionally, you can define an exact legend position using LEGEND_POSITION_X and LEGEND_POSITION_Y. To use automatic legend positioning, do NOT define these values
#Optionally, shape scaling can be present (LEGEND_SHAPE_SCALES). For each shape, you can define a scaling factor between 0 and 1.
#To order legend entries horizontally instead of vertically, set LEGEND_HORIZONTAL to 1
#Shape should be a number between 1 and 6, or any protein domain shape definition.
#1: square
#2: circle
#3: star
#4: right pointing triangle
#5: left pointing triangle
#6: checkmark

#LEGEND_TITLE Dataset_legend
#LEGEND_SCALE 1
#LEGEND_POSITION_X 100
#LEGEND_POSITION_Y 100
#LEGEND_HORIZONTAL 0
#LEGEND_SHAPES 1 1 2 2
#LEGEND_COLORS #ff0000 #00ff00 rgba(0,255,0,0.5) #0000ff
#LEGEND_LABELS value1 value2 value3 value4
#LEGEND_SHAPE_SCALES 1 1 0.5 1

#width of the colored strip
#STRIP_WIDTH 25

#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap.
#MARGIN 0

#border width; if set above 0, a border of specified width (in pixels) will be drawn around the color strip 
#BORDER_WIDTH 0

#border color; used when BORDER_WIDTH is above 0
#BORDER_COLOR #0000ff

#if set to 1, border will be drawn completely around each colored strip box
#COMPLETE_BORDER 0

#always show internal values; if set, values associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
#SHOW_INTERNAL 0


#display or hide the individual label inside each colored strip (when defined in the data below)
#SHOW_STRIP_LABELS 1

#position of the strip label within the box; 'top', 'center' or 'bottom'
#STRIP_LABEL_POSITION center

#strip label size factor (relative to the tree leaf labels)
#STRIP_LABEL_SIZE_FACTOR 1


#rotation of the strip labels; used only in rectangular tree display mode
#STRIP_LABEL_ROTATION 0

#strip label shift in pixels (positive or negative)
#STRIP_LABEL_SHIFT 0

#STRIP_LABEL_COLOR #000000

#draw a black outline around the text (width in pixels)
#STRIP_LABEL_OUTLINE 0.5

#calculate the label color automatically (black or white), based on the darkness of the color strip
#STRIP_LABEL_AUTO_COLOR 0

#display or hide the dataset label above the colored strip
#SHOW_LABELS 1

#dataset label size factor
#SIZE_FACTOR 1

#dataset label rotation
#LABEL_ROTATION 0

#dataset label shift in pixels (positive or negative)
#LABEL_SHIFT 0

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages

#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA

#Examples:
#assign a red colored strip to leaf 9606, with label 'Human'
#9606 #ff0000 Human

#assign a green, semi-transparent (alpha 0.5) strip to an internal node, without any label. If 'Show internal values' is set to 'No', this will only be displayed if the node is collapsed. 
#9606|5664 rgba(0,255,0,0.5)\n"""

    def generate_random_color_code(self):
        # Generate three random numbers between 0 and 255
        R = random.randint(0, 255)
        G = random.randint(0, 255)
        B = random.randint(0, 255)
        
        # Convert the RGB values to a hexadecimal color code
        color_code = "#{:02X}{:02X}{:02X}".format(R, G, B)
        
        return color_code
    
    def write_fo(self, df, colors):
        # Output file format
        # Leaf_ID   Leaf_Color  Leaf_Label
        # 123   #456727 number
        # hello #C23B57 str
        raw_df = {"Leaf_ID": df["Leaf_ID"][:],
                  "Leaf_Color": [],
                  "Leaf_Label": df["Leaf_Label"][:]}

        for elm in raw_df["Leaf_Label"]:
            raw_df["Leaf_Color"].append(colors[elm])

        new_df = pd.DataFrame(raw_df)
        
        # Color table in csv file
        print('[*] Generating color_table.csv...')
        new_df.to_csv(os.path.join(self.fo, "color_table.csv"), sep='\t', index=False)

        # Annotation file
        print('[*] Generating color_strip.txt...')
        with open(os.path.join(self.fo, "color_strip.txt"), 'w') as out:
            out.write(f'{self.annotation_prefix}\n')
            for i in range(new_df.shape[0]):
                out.write(f'{new_df["Leaf_ID"][i]} {new_df["Leaf_Color"][i]} {new_df["Leaf_Label"][i]}\n')

    # main
    def generator(self):
        # Input file format
        # Leaf_ID   Leaf_Label
        # 123   number
        # hello str
        print('[*] Reading data...')
        df = pd.read_csv(self.fi, sep=self.sep)
        
        # Labels
        labels = set()
        for elm in df["Leaf_Label"]:
            labels.add(elm)
        
        # Add color
        print('[*] Assigning color...')
        colors = {}
        color_used = []
        for elm in labels:
            if elm.lower()  == 'fibroblast':
                colors[elm] = '#f44336'
            
            else:
                while True:
                    color = self.generate_random_color_code()

                    if color not in color_used:
                        colors[elm] = color
                        break
        
        print("[*] Generating output files...")
        self.write_fo(df, colors)
        print('[*] Done.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-I', '--FILE_INPUT', required=True, type=str)
    parser.add_argument('-O', '--FILE_OUTPUT_PATH', required=True, type=str, help="File name is not required.")
    parser.add_argument('-SEP', '--SEPARATOR', default="\t", help="Default is TAB.")
    args = parser.parse_args()
    fi = args.FILE_INPUT
    fo = args.FILE_OUTPUT_PATH
    sep = args.SEPARATOR

    annotate = iTOL_annotation(fi, fo, sep)
    annotate.generator()
