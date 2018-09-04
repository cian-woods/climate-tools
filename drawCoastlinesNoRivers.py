def drawCoastlinesNoRivers(m,color='k',linewidth=1):
    poly_stop = 85
    fill_color = '0.9'
    # If you don't want lakes set lake_color to fill_color
    #m.fillcontinents(color=fill_color,lake_color=fill_color)
    # Draw the coastlines, with a thin line and same color as the continent fill.
    coasts = m.drawcoastlines(zorder=1,color='white',linewidth=0)
    # Exact the paths from coasts
    coasts_paths = coasts.get_paths()
    # In order to see which paths you want to retain or discard
    # you'll need to plot them one at a time noting those that you want etc.
    for ipoly in range(0,75+1,1):
      if ipoly > poly_stop: continue
      r = coasts_paths[ipoly]
      # Convert into lon/lat vertices
      polygon_vertices = [(vertex[0],vertex[1]) for (vertex,code) in
                          r.iter_segments(simplify=False)]
      px = [polygon_vertices[i][0] for i in xrange(len(polygon_vertices))]
      py = [polygon_vertices[i][1] for i in xrange(len(polygon_vertices))]
      m.plot(px,py,'k-',linewidth=linewidth,color=color,zorder=50)
      #pl.title(ipoly)
      #pl.show()
      #pl.clf()
