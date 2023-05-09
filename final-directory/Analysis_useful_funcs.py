import numpy as np
import math
import matplotlib.pyplot as plt
import os
import matplotlib

sKy_colors = {'light blue':'#63B8FF', 'blue':'#4876FF', 'very dark blue':'#27408B', 'blue grey':'#C6E2FF', 'dim cyan':'#98F5FF', 'cyan':'#00FFFF','red':'#FF4040', 'mute red':'#EE6363', 'dark mute red':'#CD5555', 'dark red':'#CD2626', 'green':'#00FF7F', 'honest green':'#008B45', 'dark green':'#008B45', 'grey':'#8B8989', 'dark grey':'#666666', 'orange':'#FF9912', 'purple':'#8E388E', 'magenta':'#FF00FF', 'purple pink':'#FF83FA', 'dark purple pink':'#BF3EFF', 'bright brown':'#8B5A00', 'dull brown':'#8B4726', 'mute brown':'#BC8F8F'}

font1 = 'Shree Devanagari 714' #use this for unbolded text
font2 = 'hiragino sans' #use this for bolded text
matplotlib.rcParams["font.family"] = font1
matplotlib.rc('axes', unicode_minus=False)
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = font1
matplotlib.rcParams['mathtext.it'] = font1
matplotlib.rcParams['mathtext.bf'] = font1

def plot_available_fonts(save_loc = 'fonts.png', bold = False):
    import matplotlib.font_manager
    fpaths = matplotlib.font_manager.findSystemFonts()

    all_fonts = []

    for i in fpaths:
        
        try:
            f = matplotlib.font_manager.get_font(i)
            print(f.family_name)
            all_fonts.append(f.family_name)
        except:
            pass
    W = 300
    H = 20
    xs_added = []
    values_to_skip = []
    #values_to_skip = [1, 2, 4, 5, 7, 8, 9, 11, 12, 13, 14, 15, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 29, 30, 31, 32, 37, 41, 43, 44, 48, 49, 50, 52, 55, 56, 59, 62, 65, 67, ]
    plt.figure(figsize=(W, H))
    for i in range(len(all_fonts)):
        if i in values_to_skip:
            continue
        try:
            xs_added.append((math.floor(i*0.1))*160.0)
            if bold:
                plt.text(x = (math.floor(i*0.1))*160.0, y = (i*0.1)%1, s = all_fonts[i]+": Absolute Magnitude", font_properties = all_fonts[i], fontsize = 20, weight = 'bold')
            else:
                plt.text(x = (math.floor(i*0.1))*160.0, y = (i*0.1)%1, s = all_fonts[i]+": Absolute Magnitude", font_properties = all_fonts[i], fontsize = 20)
            plt.text(x = (math.floor(i*0.1))*160.0, y = (i*0.1)%1+0.03, s = str(i), fontsize = 15)
        except:
            pass

    plt.xlim([0, 10000])
    plt.savefig(save_loc)

def isint(element):
    return is_int(element)
    
def exists(file):
    return os.path.exists(file)
    
def is_int(element):
    '''checks if an element can be turned into an int'''
    try:
        int(element)
        return True
    except ValueError:
        return False

def isfloat(element):
    return is_float(element)
    
    
def conv_to_deluxe_table(file_loc, deluxe_table_loc = None, table_caption = None, table_label = None, separators = [' ', ','], insert_columns = None, switch_columns = None, insert_rows = None):
    '''inputs: file_loc-location of file with data to turn into latex code for deluxe table, deluce_table_loc-location of where to write code for deluxe table to, insert_columns-l-the number of columns to be inserted-this will create this many columns near the end of the table, with randomized column names, switch_columns-will move columns around-this should be a list like["1,2"]-which will tell it to switch columns 1 and 2-indexing from 0, insert_rows-number of rows to be inserted-this will be done after all the column switching and adding is done
       outputs: none
       this will take a location of a file that holds data and convert it into latex code for a deluxe table. If deluxe_table_loc is passed, the code will be written to that location
    '''
    
    random_col_names = ['Cheese', 'Lisa', 'Jennie', 'Pizza', 'Bibimbap', 'PaniPuri', 'PavBhaji', 'Fries', 'Dip', 'Hummus', 'Pita']
    data = read_file(file_loc, separators = separators)['dict']
    
    keys = list(data.keys())
    nrow = len(data[keys[0]])
    if insert_columns is not None:
        for i in range(insert_columns):
            if i+len(keys) < len(random_col_names)-1:
                colname = random_col_names[i+len(keys)]
            else:
                colname = "Column "+str(i)
            data[colname] = [str(i) for i in np.zeros(nrow)]
            
    keys = list(data.keys())
    
    for switch in switch_columns:
        a = keys[int(switch[:switch.index(',')])]
        b = keys[int(switch[switch.index(',')+1:])]
        
        temp = [i for i in data[a]]
        data[a] = [i for i in data[b]]
        data[b] = [i for i in temp]
    
    keys = list(data.keys())
    ncol = len(keys)
    
    latex_code = "\\begin{deluxetable}{"+"c"*ncol+"}\n"
    latex_code += "\\tabletypesize{\\footnotesize}\n"
    latex_code += "\\tablecolumns{"+str(ncol)+"}\n"
    
    if table_caption is None:
        latex_code += "\\tablecaption{I like Cheese. "
    else:
        latex_code += "\\tablecaption{"+str(table_caption)+". "
    if table_label is None:
        latex_code += "\\label{table:fun_table}.}\n"
    else:
        latex_code += "\\label{table:"+str(table_label)+"}.}\n"
        
    col_heads = ""
    for i in range(ncol-1):
        if i > len(random_col_names)-1:
            col_name = keys[i]
        else:
            col_name = random_col_names[i]
        col_heads += " \\colhead{"+str(col_name)+"} &"
        
    if i > len(random_col_names)-1:
        col_name = keys[i]
    else:
        col_name = random_col_names[i]
    
    col_heads += " \\colhead{"+str(col_name)+"} "
    
    latex_code += "\\tablehead{"+str(col_heads)+"}\n"
    latex_code += "\startdata\n"
    
    for i in range(nrow):
        #read the value at each row from the keys
        line = ''
        for j in range(len(keys)-1):
            key = keys[j]
            line += data[key][i] + " & "
        #print("keys: "+str(keys))
        #print("keys[len(keys)-1]: "+str(keys[len(keys)-1]))
        line += data[keys[-1]][i]+" \\\\\n"
        latex_code += line
        
    latex_code += '\enddata\n'
    latex_code += '\end{deluxetable}'
    
    if deluxe_table_loc is not None:
        write_to_file(deluxe_table_loc, latex_code)
    else:
        print(latex_code)
    
    
    
    
    
def ls(dir = None, extension = None):
    '''inputs: dir-directory to read the list of files from, extension-which extension should the files have to be listed (default to none-will list all file types)
       outputs: list of files WITHOUT directory name
       This function gets a list of all files within a directory name. Optionally, can input a specific file extension to only list out files with that extension. This will only list out the short name of the files, like ls would on terminal
    
    '''
    
    files_list = []
    if dir is not None:
    
    
        for filename in os.listdir(dir):
            ext = filename[-1*len(extension):]
            short_name = filename
            if extension is not None and ext == extension:
                files_list.append(short_name)
    else:
    
        for filename in os.listdir():
            ext = filename[-1*len(extension):]
            short_name = filename
            if extension is not None and ext == extension:
                files_list.append(short_name)
        
    return files_list
            


def write_to_file(file_loc, data_out, separator = ' ', headers=None, append=False):
    '''inputs: file_loc-location to which to write to, data_to_write-this is a 2d array that will be written to a file
       outputs: none
       This function writes a 2d array to a file'''
       
    if not append:
        out = open(file_loc ,'w')
    else:
        out = open(file_loc ,'a+')
    
    if type(data_out) is str:
        out.write(data_out+"\n")
        
    else:
        for i  in range(len(data_out)):
            if isinstance(data_out[i], float) or type(data_out[i]) is str:
                out.write(str(data_out[i])+separator)
            else:
                for j in range(len(data_out[i])):
                    out.write(str(data_out[i][j])+separator)
                
            out.write("\n")
    out.close()
            
            

def load_file(file_loc, header_marker=None):
    '''same as read_file'''
    
    return read_file(file_loc, header_marker = header_marker)

def read_file(file_loc, comments ='#', separators = [], skip_lines = [-1], header_list = None, header_row = None, header_line = None, header_marker=None, find_floats = False):
    '''inputs: file_loc-string of file location to read in, comments-which character makes the line a comment, separators-list of characters that are used to separate line entry from line entry, skip_lines-list of lines to skip when reading in, header_row-which row in the datafile is the set of headers (defaults to none, if none, the dictionary will have keys of just row numbers), header_line-aame as header_row, header_marker-the character that marks a line as the header, find_floats-will only return a data entry in the final dictionary if it is a float
       outputs: data-2d array of data read in and passed out
       This function reads in a dataset in a file located at file_loc and returns a dictionary of this'''
       
    if header_line != None and header_row == None:
        header_row = header_line
        
        
    
        
    headers = ''
    if header_list is not None:
        headers = header_list.rstrip()
    row_num = 0
    it = 0
    skipped_lines = []
    with open(file_loc) as f:
        lines = []
        for line in f:
            if row_num not in skip_lines:
                if len(separators) > 0:
                    for separator in separators:
                        line = line.replace(separator, ' ')
                if len(line.split()) > 0 and line.split()[0] != comments:
                    if len(headers) == 0 and header_row is not None and it == 0:
                        headers = line.rstrip()
                        
                    else:
                        lines.append(line.rstrip().split())
                    it += 1
            elif row_num in skip_lines:
                skipped_lines.append(line)
            row_num += 1
        
    if len(separators) > 0:
        for separator in separators:
            headers = headers.replace(separator, ' ')
    
    if len(headers) > 0:
        headers = headers.split()
    
    dict = {}
    if len(headers) > 0:
        for header in headers:
            dict[header] = []
    
    for i in range(len(lines)):
        for j in range(len(lines[i])):
            if len(headers) < j+1:
                if j not in dict.keys():
                    dict[j] = []
                dict[j].append(lines[i][j])
            else:
                if headers[j] not in dict.keys():
                    dict[headers[j]] == []
                dict[headers[j]].append(lines[i][j])
        
    return {"lines": lines, "dict": dict, "skipped_lines":skipped_lines}
        
def correct_reddening(dataset, wavelength_angstroms, EBV, band = None, assume_optical = True):
    '''corrects a given band of observations given EBV extinction along that line of sight. This uses the equation: M_true = M_observed + A, derived from f_lambda^true = f_lambda^obs*10^(A/2.5)'''
    '''If band is given, then the standard R value for that band is used to do the extinction correction'''
    
    lambda_alpha = 6563
    lambda_beta = 4861
    R_h_alpha = 2.535
    R_h_beta = 3.609
    #I know R goes as 1/lambda in optical, so
    #R = k/lambda, using K = R1*lambda1 = R2 * lambda2, using the value at halpha, you get K = 16637.2 and at hbeta gives K = 17543.3, so I'll use the average K = 17090.3
    band_standard_R_values = {"B" : 4.1, "V" : 3.1, "R":2.3, "I":1.5}
    if band is None or band.upper() not in band_standard_R_values.keys():
        K = 17090.3
        R_optical = 17090.3/wavelength_angstroms
        
        if not(wavelength_angstroms > 2000 and wavelength_angstroms < 9000):
            assume_optical = False
        if assume_optical:
            R = R_optical
            
    else:
        R = band_standard_R_values[band.upper()]
        
    A = R*EBV
    
    #notice: R_B is 4.1 versus R_V is 3.1. We also know that dust extinction causes 'reddening'. As in, the flux in B band is decreased more than the flux in V band. Therefore, the B-mag value is numerically increased more than the V-mag value. Therefore, A is a number that is added to each band's magnitude by the dust. To correct for the dust, we must subtract A away from each band.
    if isfloat(dataset):
        return dataset - A
    else:
        return np.asarray(dataset)-A


def is_float(element):
    '''checks if an element can be turned into a float'''
    try:
        float(element)
        return True
    except:
        return False


def find_nearest(wv,value):
    '''inputs: wv-list of data to look in, value-the value we are looking for in wv
       outputs: idx-index in wv that is closest to value
       this finds the index of the closest element to value in list wv'''
    idx= (np.abs(np.asarray(wv) - value)).argmin()
    return idx

def mkdir(dir_loc):
    '''inputs: dir_loc-location of directory to create
       outputs: none
       This function creates the directory at dir_loc if it doesn't exist. If it does exist, nothing happens'''
    if not os.path.isdir(dir_loc):
        os.mkdir(dir_loc)
        
def get_order_of_mag(num):
    '''returns order of magnitude for a number passed in. This returns just the exponent, so if you pass in 8.9E-45, it will return -45'''
    
    if num == 0:
        return 0

    absnum = abs(num)
    order = math.log10(absnum)
    res = math.floor(order)

    return res
    

def pretty_plot(x, y, xlabel = r'', ylabel = '', title = '', label = '', xlim = [], ylim = [], xscale='', yscale='', save_loc='', display_or_nah = False, scatter = False, s = None, dark_theme = False, minimalist = False, color = None, bold = False):
    '''inputs: x-indpendent variable of plot to create, y-dependent variable of plot to create, xlabel-string of label on x-axis, ylabel-string of label on y-axis, title-string of title for figure, xlim-list of limits of plot for x-axis, ylim-list of limits of plot for y-axis, labels-string of label of plot, save_loc-string location in file where to store the plot (default to ''; if '' is set, the figure won't be saved), display_or_nah-boolean variable to control whether the figure should be displayed (default to False)
        outputs: none
        This function plots x and y using plt.plot--I only have this function because I have a specific way I like plotting, and hate remembering all the details, so this function plots things for me'''
    
    
    weight = 'normal'
    
    if bold:
        plt.rcParams["font.family"] = font2
        weight = 'bold'
    else:
        plt.rcParams["font.family"] = font1
    plt.rc('axes', unicode_minus=False)
    W = 15
    H = 10

    plt.figure(figsize=(W, H))
    
    '''
    if dark_theme:
        plt.style.use('dark_background')
    elif minimalist:
        plt.style.use('fivethirtyeight')
    else:
        plt.style.use('seaborn')'''
        
    
    if color is None:
        color = sKy_colors[list(sKy_colors.keys())[0]]
    
    if not scatter:
    
        if len(label) > 0:
            plt.plot(x, y, linewidth=H/3, label=label, color = color)
        else:
            plt.plot(x, y, linewidth=H/3, color = color)
    else:
        if s == None:
            s = 5
        if len(label) > 0:
            plt.scatter(x, y, s = H*s, label=label, color = color)
        else:
            plt.scatter(x, y, s = H*s, color = color)
    if len(xlabel) > 0:
        plt.xlabel(xlabel, fontsize=2*W, weight = weight)
        
    plt.xticks(fontsize=1.5*W)
    
    if len(xlim) > 0:
        plt.xlim(xlim)
        
    if len(ylim) > 0:
        plt.ylim(ylim)
    
    if len(ylabel) > 0:
        plt.ylabel(ylabel, fontsize=2*W, weight = weight)
        
    plt.yticks(fontsize=1.5*W)
    
    if len(title) > 0:
        plt.title(title, fontsize=2*W, weight = weight)
        
    if len(label):
        plt.legend(fontsize = 1.5*H, weight = weight)
    
    if xscale == 'log':
        plt.xscale('log')
        
    if yscale == 'log':
        plt.yscale('log')
        
    if len(save_loc) > 0:
        plt.savefig(save_loc)
        
    if display_or_nah:
        plt.show()
        
    plt.close()
        
    return plt
