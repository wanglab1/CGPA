###################################################################################
################ Correlated genes tab functions ########################
#####################################################################################
corr_ui = function(input_gene,cancer_type,cancer_options,action_button,data_output,plot_output){
  #-------------------------------------------------------------------------
  # This func is for correlated mRNA and lncRNA correlations 
  # input_gene
  # cancer_type
  # cancer_options
  # action_button
  # data_output
  # plot_output
  #-----------------------------------------------------------------------
           fluidRow(
             br(),
             column(2,  
                    sel_input_ui(cancer_type,"Cancer type", choices = cancer_options,selected=c("LUSC")),
                    actionButton(action_button, "Run")
                    
             ),
             br(),
             column(5,
                    div(id="box1",
                        box(
                          title="Top mRNA corrlated with the selected lncRNA (Table)",collapsible = TRUE,status="warning",
                          solidHeader = TRUE,
                          width = 12,height=600,
                          withLoader(DTOutput(data_output),type="html",loader="loader1")
                        )
                    )
             ),
             column(5,
                    div(id="box1",
                        box(
                          title="Top mRNA corrlated with the selected lncRNA (Figure)",collapsible = TRUE,status="warning",
                          solidHeader = TRUE,
                          width = 12,height=600,
                          withLoader(plotOutput(plot_output,height="600px"),type="html",loader="loader1")
                        )
                    )
             )
           )
  
}


