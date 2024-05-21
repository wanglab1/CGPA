#=================================================================================
####################### Pan cancer dashboard tab ################################
###################### 06/29/2022 ##############################
#=================================================================================
pan_ui = function(box_id,title,height,output_fig){ 
#-------------------------------------------------------------
# box_id: box_id for each panel, different box has diff styles
# title: titles for each box
# height: box height
# output figure name
#-------------------------------------------------------------
  div(id=box_id,
      box(
        title=title,collapsible = TRUE, 
        status="warning",
        solidHeader = TRUE,
        width = 12,height=height,

        if(box_id=="box4"){
          div(style='max-width: 100%; height: 500px;',
              div(style="text-align: left;",uiOutput(output_fig))

          )
        }else{
            (div(style='width: auto;max-width: 100%; height: 500px; ',
                 div(style="text-align: center;width: auto;",imageOutput(output_fig)),
            ))
          
        }

      )
  )
}

pan_server_fig = function(box_id,file=NULL,file_name=NULL,input_gene){
#----------------------------------------------------------------------
  # functions to import images
  # box_id: to identify functions in different panel
  # file: file path
  # file_name: file names
  # input_gene: geneNames
#----------------------------------------------------------------------
  if(box_id=="box4"){
    src=paste0("https://string-db.org/api/svg/network?identifiers=%0D",input_gene,"&limit=10&network_flavor=evidence&species=9606")
    src
  }else{
    tryCatch({
      file = file
      fig = paste0(file,make.names(input_gene),file_name)
      list(src=fig,alt="fig")
    },error=function(e){
      ""
    })
  }

}


pan_server_text = function(input_data,type){
  #----------------------------------------------------------------------
  # functions to summarize forest plot
  # input_data: dataset for the selected gene (from pan_server_text_data())
  # type: OS or PFI
  #----------------------------------------------------------------------
  
  uni_cox1 = input_data
  names_sig = uni_cox1$names[which(as.numeric(uni_cox1$p.value) < 0.05)]
  
  if(length(names_sig)==0){
    text= "The prognostic marker is not significant in any cancer types"
  } else if(length(names_sig)!=0){
    names_fav = uni_cox1$names[which(uni_cox1$symbol=="<1"&uni_cox1$p.value<0.05)]
    names_unfav = uni_cox1$names[which(uni_cox1$symbol==">=1"&uni_cox1$p.value<0.05)]
    if(length(names_unfav)==0){
      text= paste0("The prognostic marker is significant and favourable in ",paste(names_fav,collapse=" ")," (hazard ratio < 1, low risk of death)")
    }else if(length(names_fav)==0){
      text= paste0("The prognostic marker is significant and unfavourable in ",paste(names_unfav,collapse=" ")," (hazard ratio > 1, high risk of death)")
    }else{
      text= paste0("The prognostic marker is significant and favourable in ",paste(names_fav,collapse=" ")," (hazard ratio < 1, low risk of death)",
                   ", and unfavourable in ",paste(names_unfav,collapse=" ")," (hazard ratio > 1, high risk of death)"," for ",type)
      
    }
  }else{
    text = "Too many zeros in the selected gene"
  }
  text
}

