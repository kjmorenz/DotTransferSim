def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 50, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    num = 0
    if iteration > total:
        num = int(iteration/total)
        iteration = iteration%total
        #suffix = suffix + " + " + str(num*100) + "\%"
    percentval = 100 * (iteration / float(total))+num*100
    percent = ("{0:." + str(decimals) + "f}").format(percentval)
    filledLength = int(length * iteration // total)
    if percentval % 100 ==0:
        filledLength = length
    bar = fill * filledLength + '-' * (length - filledLength)
    
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()