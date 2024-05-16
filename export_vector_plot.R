export_vector_plot = function(plist, fileout, bg = "white", pointsize = 12, editable = TRUE, fonts = "arial") {
    doc <- read_pptx()
    for(pp in 1 : length(plist)) {
        myplot2 =  dml(ggobj = plist[[pp]],
                       bg = bg,
                       pointsize = pointsize,
                       editable = editable, 
                       fonts = fonts)
        
        
        doc <- add_slide(doc, "Title and Content", "Office Theme")
        doc <- ph_with(doc, myplot2, location = ph_location_fullsize())
    }
    
    print(doc, target = fileout)
    
}
