import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy
from scipy import io
import glob

def collector(path):
    '''Collects all IDL save data from a given path and stores each file as an element in a list.'''
    
    filenames = glob.glob(path)
    
    data = {'data':[scipy.io.readsav(filenames[i],python_dict=True) \
            for i in range(len(filenames))],'filenames':filenames}
    return data

def separator(data):
    '''Compiles data into separate lists of extended and point sources'''

    point_data = [[data['data'][i]['source_array'][j] \
        for j in range(len(data['data'][i]['source_array'])) \
        if data['data'][i]['source_array'][j][-2] is None ] \
    for i in range(len(data['data']))]
    
    extended_data = [[data['data'][i]['source_array'][j] \
        for j in range(len(data['data'][i]['source_array'])) \
        if data['data'][i]['source_array'][j][-2] is not None ] \
    for i in range(len(data['data'])) ]

    return {'extsources':extended_data,'psources':point_data}

def pixelate(ra_zoom, dec_zoom, n_bins, ra_total, dec_total, flux_total):
    import numpy as np
    #Check to see which dimension is larger so that a square in ra,dec can 
    #be returned
    if (ra_zoom[1]-ra_zoom[0]) > (dec_zoom[1]-dec_zoom[0]):
        zoom = ra_zoom
    else:
        zoom = dec_zoom

    #Find the size of the bins using the largest dimension and the num of bins
    binsize = (zoom[1]-zoom[0])/n_bins

    #Create arrays for ra and dec that give the left side of each pixel
    ra_bin_array = (np.array(range(n_bins)) * binsize) + ra_zoom[0]
    dec_bin_array = (np.array(range(n_bins)) * binsize) + dec_zoom[0]
    #Create an empty array of pixels to be filled in the for loops
    pixels = np.zeros((len(ra_bin_array),len(dec_bin_array)))

    #Histogram components into ra bins
    ra_histogram = np.digitize(ra_total,ra_bin_array)
    ###print ra_histogram

    #Begin for loop over both dimensions of pixels, starting with ra
    for bin_i in range(len(ra_bin_array) - 2):
        ###print range(len(ra_bin_array) -2
        ###print "bin_i",bin_i
        #Find the indices that fall into the current ra bin slice
        ra_inds = np.where(ra_histogram == bin_i)
        ###print "rainds", ra_inds[0]
        ###print "lenrainds", len(ra_inds[0])

        #Go to next for cycle if no indices fall into current ra bin slice
        if len(ra_inds[0]) == 0:
            continue

        #Histogram components that fall into the current ra bin slice by dec
        #print "dectotindex", dec_total[ra_inds]
        #print "decbin", dec_bin_array
        dec_histogram = np.digitize(dec_total[ra_inds],dec_bin_array)
        #print "dechist",dec_histogram
        #Begin for loop by dec over ra bin slice
        for bin_j in range(len(dec_bin_array) -2):
            
            #Find the indicies that fall into the current dec bin
            dec_inds = np.where(dec_histogram == bin_j)

            #Go to next for cycle if no indices fall into current dec bin			
            if len(dec_inds[0]) == 0:
                continue
            #Sum the flux components that fall into current ra/dec bin
            ###print "bi",bin_i,bin_j
            ###print "inds",ra_inds, dec_inds
            pixels[bin_i,bin_j] = np.sum(flux_total[ra_inds[0][dec_inds][0]])

    #Find the pixel centers in ra/dec for plotting purposes
    ra_pixel_centers = (np.arange(n_bins) * binsize) + ra_zoom[0] + binsize/2.
    dec_pixel_centers = (np.arange(n_bins) * binsize) + dec_zoom[0] + binsize/2.

    return pixels, ra_pixel_centers, dec_pixel_centers

def plotEO(data,minI,sumI,n_bins):
    from matplotlib.colors import LogNorm
    import matplotlib.patches as patches
    
    separated = separator(data)
    
    indexed_point_sources_RA = [[separated['psources'][i][j]['RA'] \
        for j in range(len(separated['psources'][i]))] \
            for i in range(len(separated['psources'])) ]

    indexed_point_sources_DEC = [[separated['psources'][i][j]['DEC'] \
        for j in range(len(separated['psources'][i]))] \
            for i in range(len(separated['psources'])) ]

    indexed_point_sources_I = [[separated['psources'][i][j]['FLUX']['I'][0] \
        for j in range(len(separated['psources'][i]))] \
            for i in range(len(separated['psources'])) ]

    indexed_EO_sources_RA = [[[separated['extsources'][i][j]['EXTEND']['RA'][k] \
        for k in range(len(separated['extsources'][i][j]['EXTEND']['RA']))] \
            for j in range(len(separated['extsources'][i]))] \
                for i in range(len(separated['extsources'])) ]

    indexed_EO_sources_DEC = [[[separated['extsources'][i][j]['EXTEND']['DEC'][k] \
        for k in range(len(separated['extsources'][i][j]['EXTEND']['DEC']))] \
            for j in range(len(separated['extsources'][i]))] \
                for i in range(len(separated['extsources'])) ]

    indexed_EO_sources_I = [[[separated['extsources'][i][j]['EXTEND']['FLUX'][k]['I'][0] \
        for k in range(len(separated['extsources'][i][j]['EXTEND']['FLUX'])) ]
            for j in range(len(separated['extsources'][i]))] \
                for i in range(len(separated['extsources'])) ]
    
    all_RA = [[indexed_point_sources_RA[i][j] \
        for j in range(len(indexed_point_sources_RA[i]))]
        + [indexed_EO_sources_RA[i][j][k] \
        for j in range(len(indexed_EO_sources_RA[i])) \
        for k in range(len(indexed_EO_sources_RA[i][j]))]
            for i in range(len(data['data'])) ]
    
    all_DEC = [[indexed_point_sources_DEC[i][j] \
        for j in range(len(indexed_point_sources_DEC[i]))]
        + [indexed_EO_sources_DEC[i][j][k] \
        for j in range(len(indexed_EO_sources_DEC[i])) \
        for k in range(len(indexed_EO_sources_DEC[i][j]))]
            for i in range(len(data['data'])) ]
    
    all_I = [[indexed_point_sources_I[i][j] \
        for j in range(len(indexed_point_sources_I[i]))]
        + [indexed_EO_sources_I[i][j][k] \
        for j in range(len(indexed_EO_sources_I[i])) \
        for k in range(len(indexed_EO_sources_I[i][j]))]
            for i in range(len(data['data'])) ]
    
    #Correcting RA to go from -180 degrees to 180 degrees instead of 0 degrees to 360 degrees.
    for i in range(len(all_RA)):
        for j in range(len(all_RA[i])):
            if all_RA[i][j] > 180:
                all_RA[i][j] -= 360    

    
    for i in range(len(data['data'])):
        
        semi_zoom = 3
        pixelreplacement = 1e-3
        
        all_ra_zoom = [min(all_RA[i]),max(all_RA[i])]
        all_dec_zoom = [min(all_DEC[i]),max(all_DEC[i])]
        all_n_bins = n_bins
        all_ra_total = np.array(all_RA[i])
        all_dec_total = np.array(all_DEC[i])
        all_flux_total = np.array(all_I[i])
        
        (all_pixels, all_ra_pixel_centers, all_dec_pixel_centers) = \
                pixelate(all_ra_zoom,all_dec_zoom,all_n_bins,all_ra_total,all_dec_total,all_flux_total)
        
        all_pixels[all_pixels == 0] = pixelreplacement
        all_logpixels = np.log10(all_pixels)
        
        for j in range(len(separated['extsources'][i])):
            if (max(indexed_EO_sources_I[i][j]) > minI) or (sum(indexed_EO_sources_I[i][j]) > sumI):
                
                
                
                EO_ra_zoom = [min(indexed_EO_sources_RA[i][j]),max(indexed_EO_sources_RA[i][j])]
                EO_dec_zoom = [min(indexed_EO_sources_DEC[i][j]),max(indexed_EO_sources_DEC[i][j])]
                EO_n_bins = n_bins
                EO_ra_total = np.array(indexed_EO_sources_RA[i][j])
                EO_dec_total = np.array(indexed_EO_sources_DEC[i][j])
                EO_flux_total = np.array(indexed_EO_sources_I[i][j])
                
                (EO_pixels, EO_ra_pixel_centers, EO_dec_pixel_centers) = \
                pixelate(EO_ra_zoom,EO_dec_zoom,EO_n_bins,EO_ra_total,EO_dec_total,EO_flux_total)
                
                EO_pixels[EO_pixels == 0] = pixelreplacement
                EO_logpixels = np.log10(EO_pixels)
                
                semi_ra_zoom = [min(indexed_EO_sources_RA[i][j])-semi_zoom,max(indexed_EO_sources_RA[i][j])+semi_zoom]
                semi_dec_zoom = [min(indexed_EO_sources_DEC[i][j])-semi_zoom,max(indexed_EO_sources_DEC[i][j])+semi_zoom]
                semi_n_bins = n_bins
                semi_ra_total = np.array(all_RA[i])
                semi_dec_total = np.array(all_DEC[i])
                semi_flux_total = np.array(all_I[i])

                (semi_pixels, semi_ra_pixel_centers, semi_dec_pixel_centers) = \
                        pixelate(semi_ra_zoom,semi_dec_zoom,semi_n_bins,semi_ra_total,semi_dec_total,semi_flux_total)

                semi_pixels[semi_pixels == 0] = pixelreplacement
                semi_logpixels = np.log10(semi_pixels)
            
                #cmap = matplotlib.cm.get_cmap('afmhot')
                cmap = matplotlib.cm.get_cmap('gist_heat')

                #cmap.set_bad((0,0,0))
                
                ############################################################################

                fig = plt.figure(figsize=(18,4))
                fig.suptitle('ObsID: {} at Frequency {} MHz'.format\
                    ([int(s) for s in re.findall('\d+',data['filenames'][i])][0],\
                    separated['extsources'][i][j]['FREQ']), fontsize = 20)

                ax = fig.add_subplot(1,3,1)
                
                all_plot = ax.imshow(np.transpose(all_logpixels), \
                    origin = "lower", interpolation = "nearest", cmap = cmap,\
                    extent = [all_ra_pixel_centers[0], all_ra_pixel_centers[-1], \
                    all_dec_pixel_centers[0], all_dec_pixel_centers[-1]])
                ax.add_patch(patches.Rectangle((min(EO_ra_pixel_centers),min(EO_dec_pixel_centers)), \
                    semi_zoom,semi_zoom,fill=False,color='cyan'))
                ax.set_xlabel('RA', fontsize = 12)
                ax.set_ylabel('DEC', fontsize = 12)
                ax.tick_params(size = 6, labelsize = 10)
                ax.minorticks_on()
                ax.tick_params('both', length=8, width=1.8, which='major')
                ax.tick_params('both',length=3, width=1.4, which='minor')
                cb = fig.colorbar(all_plot)
                cb.set_label('Log Janskies')
                
                ax = fig.add_subplot(1,3,2)
                semi_plot = ax.imshow(np.transpose(semi_logpixels), \
                    origin = "lower", interpolation = "nearest", cmap = cmap,\
                    extent = [semi_ra_pixel_centers[0], semi_ra_pixel_centers[-1], \
                    semi_dec_pixel_centers[0], semi_dec_pixel_centers[-1]])
                ax.add_patch(patches.Rectangle((min(EO_ra_pixel_centers),min(EO_dec_pixel_centers)), \
                    (max(EO_ra_pixel_centers)-min(EO_ra_pixel_centers)),\
                    (max(EO_dec_pixel_centers)-min(EO_dec_pixel_centers)),fill=False,color='yellow'))
                ax.set_xlabel('RA', fontsize = 12)
                ax.set_ylabel('DEC', fontsize = 12)
                ax.tick_params(size = 6, labelsize = 10)
                ax.minorticks_on()
                ax.tick_params('both', length=8, width=1.8, which='major')
                ax.tick_params('both',length=3, width=1.4, which='minor')
                cb = fig.colorbar(semi_plot)
                cb.set_label('Log Janskies')
                
                
                ax = fig.add_subplot(1,3,3)
                EO_plot = ax.imshow(np.transpose(EO_logpixels), \
                    origin = "lower", interpolation = "nearest", cmap = cmap,\
                    extent = [EO_ra_pixel_centers[0], EO_ra_pixel_centers[-1], \
                    EO_dec_pixel_centers[0], EO_dec_pixel_centers[-1]])

                ax.set_xlabel('RA', fontsize = 12)
                ax.set_ylabel('DEC', fontsize = 12)
                ax.tick_params(size = 6, labelsize = 10)
                ax.minorticks_on()
                ax.tick_params('both', length=8, width=1.8, which='major')
                ax.tick_params('both',length=3, width=1.4, which='minor')
                ax.annotate("Source ID: {}\nComponents: {}\nHighest Flux: {}\nTotal Flux: {}\nMean RA: {}\nMean Dec: {}".format \
                    (separated['extsources'][i][j]['ID'], \
                    len(indexed_EO_sources_I[i][j]),\
                    np.max(indexed_EO_sources_I[i][j]), \
                    np.sum(indexed_EO_sources_I[i][j]), \
                    np.mean(indexed_EO_sources_RA[i][j]),\
                    np.mean(indexed_EO_sources_DEC[i][j])),\
                    xy=(2, .33), xytext=(0, 0), xycoords=('axes fraction', 'figure fraction'),\
                    textcoords='offset points', size=12, ha='center', va='bottom',\
                    bbox=dict(boxstyle="square", fc="w"))
                
                cb = fig.colorbar(EO_plot)
                cb.set_label('Log Janskies')
                #plt.savefig('pixelatedEO'+'{}'.format(separated['extsources'][i][j]['ID'])+'.png')
    return plt.show()