rm(list=ls())

#### generate data set #####
simu_data<-function(iseed,t,n,p1,tau,J,gamma){
  set.seed(iseed)
  data_example<-list()
  t<-t
  n<-n ## unique number of ID
  p1<-p1 ## treatment assignment
  #x1<-rnorm(n,1,1)
  #x2<-rnorm(n,0,2)
  tau<-tau
  J<-J ## upper for recurrent event number each person
  omega<-rgamma(1,1/tau,1/tau)
  group<-as.numeric(rbinom(n,1,p1))
  gamma<-rep(gamma,n)
  #C<-round(as.numeric(100*rexp(n,exp((x1+x2)))))+1
  #C<-round(as.numeric(rnorm(n,200,10)))
  C<-rep(365,n)
  id<-seq(1,n,by=1)
  lambda0<-0.9  ## assuming one piece constant

  for (i in 1:n){
    j<-1
    r_new<-0
    gap_time<-round(100*rexp(1,lambda0*omega*exp(group[i]*gamma[i])))+1
    s2<-ifelse(gap_time>=C[i]|gap_time>=t,0,1)
    repeat{
      if (j>=2)
      {r_new<-round(100*rexp(1,lambda0*omega*exp(group[i]*gamma[i])))+1}
      if (((sum(gap_time)+r_new)<t | (sum(gap_time)+r_new)<C[i]) & length(gap_time)<J) {gap_time=c(gap_time,r_new)
      s2=c(s2,1)}
      else if (((sum(gap_time)+r_new)>C[i] |(sum(gap_time)+r_new)>t)  & length(gap_time)<=J+1){
        gap_time<-c(gap_time,abs(min(t,C[i])-sum(gap_time)))
        #gap_time=gap_time[-1]
        #s2<-s2[-1]
        s2<-c(s2,0)
        break
      }
      else if (((sum(gap_time)+r_new)>C[i]|(sum(gap_time)+r_new)>t) & (sum(gap_time)<min(C[i],t))){
        gap_time<-c(gap_time,min(t,C[i])-sum(gap_time))
        #gap_time=gap_time[-1]
        #s2<-s2[-1]
        s2<-c(s2,0)
        break
      }
      else if(length(gap_time)>=J+1){
        gap_time=gap_time[-1]
        s2=s2[-1]
        break
      }
      j=j+1
    }




    data_example[[i]]<-cbind(id=id[i],gap_time=gap_time,group=group[i],lambda0=lambda0,c=C[i],
                             gamma=gamma[i],event=s2,t=t,days=cumsum(gap_time))
  }
  data_example<-data.frame(do.call("rbind",data_example))
  data_example<-data_example[data_example$gap_time!=0,]
  data_example$days<-ifelse(data_example$days>data_example$c,
                            data_example$c,data_example$days)

  data_example$event<-ifelse(data_example$days<data_example$c|data_example$days<data_example$t,
                             1,0)
  data_example$event<-ifelse(data_example$days>=data_example$c|data_example$days>=data_example$t,
                             0,1)
  return(data_example)
}

data_example1<-simu_data(1234,365,40,0.5,0.001,30,0)
data_example2<-simu_data(1234,365,60,0.5,0.001,30,0)
data_example3<-simu_data(1234,365,80,0.5,0.001,30,0)

data_example4<-simu_data(1234,365,40,0.5,0.001,30,-1.5)
data_example5<-simu_data(1234,365,60,0.5,0.001,30,-1.5)
data_example6<-simu_data(1234,365,80,0.5,0.001,30,-1.5)

data_example7<-simu_data(100,365,40,0.5,2,30,0)
data_example8<-simu_data(100,365,60,0.5,2,30,0)
data_example9<-simu_data(100,365,80,0.5,2,30,0)

data_example10<-simu_data(100,365,40,0.5,2,30,-1.5)
data_example11<-simu_data(100,365,60,0.5,2,30,-1.5)
data_example12<-simu_data(100,365,80,0.5,2,30,-1.5)

data_example13<-simu_data(1000,365,40,0.5,1,30,0)
data_example14<-simu_data(1000,365,60,0.5,1,30,0)
data_example15<-simu_data(1000,365,80,0.5,1,30,0)

data_example16<-simu_data(1000,365,40,0.5,1,30,-1.5)
data_example17<-simu_data(1000,365,60,0.5,1,30,-1.5)
data_example18<-simu_data(1000,365,80,0.5,1,30,-1.5)

library(reda)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("ggpubr")
library(ggpubr)
#install.packages("reReg")
library(reReg)
#install.packages("tidyverse")
library(tidyverse)
library(shiny)
library(DBI)

data_example1$term_events<-0
data_example1<-as.data.frame(data_example1)

data_example2$term_events<-0
data_example2<-as.data.frame(data_example2)

data_example3$term_events<-0
data_example3<-as.data.frame(data_example3)

data_example4$term_events<-0
data_example4<-as.data.frame(data_example4)

data_example5$term_events<-0
data_example5<-as.data.frame(data_example5)

data_example6$term_events<-0
data_example6<-as.data.frame(data_example6)

data_example7$term_events<-0
data_example7<-as.data.frame(data_example7)

data_example8$term_events<-0
data_example8<-as.data.frame(data_example8)

data_example9$term_events<-0
data_example9<-as.data.frame(data_example9)

data_example10$term_events<-0
data_example10<-as.data.frame(data_example10)

data_example11$term_events<-0
data_example11<-as.data.frame(data_example11)

data_example12$term_events<-0
data_example12<-as.data.frame(data_example12)

data_example13$term_events<-0
data_example13<-as.data.frame(data_example13)

data_example14$term_events<-0
data_example14<-as.data.frame(data_example14)

data_example15$term_events<-0
data_example15<-as.data.frame(data_example15)

data_example16$term_events<-0
data_example16<-as.data.frame(data_example16)

data_example17$term_events<-0
data_example17<-as.data.frame(data_example17)

data_example18$term_events<-0
data_example18<-as.data.frame(data_example18)


prepare_descriptive_table <- function(df, digits = c(0, 3, 3, 3, 3, 3, 3, 3), format = "html") {
  #if(!is.data.frame(df)) stop("df needs to be a dataframe")
  df <- as.data.frame(df)
  #df <- df[sapply(df, is.logical) | sapply(df, is.numeric)]
  df <- df[sapply(df, is.logical) | sapply(df, is.numeric)|sapply(df, is.character)]
  if ((ncol(df) < 1) | (nrow(df) < 2)) stop("insuitable data frame (does it contain numerical data?)")
  if (!is.numeric(digits) | length(digits) != 8) stop("digits vector is not numeric or has wrong length (!= 8)")

  t <- cbind(apply(df,2,function(x) as.integer(sum(!is.na(x)))),
             apply(df,2,mean, na.rm=TRUE),
             apply(df,2,stats::sd, na.rm=TRUE),
             t(apply(df, 2, function(x) stats::quantile(x, probs=c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))))
  colnames(t) <- c("N", "Mean", "Std. dev.", "Min.", "25 %", "Median", "75 %", "Max.")
  t <- as.data.frame(t)
  t$N <- as.integer(t$N)
  t <- t[which(!is.na(digits))]
  digits <- digits[!is.na(digits)]
  list(df = t, kable_ret = knitr::kable(t, format, digits,
                                        format.args = list(big.mark = ","),
                                        caption = "Descriptive Statistics"))
}


library(shiny)
#library(foreach)
#library(doParallel)
warning=F


data_example1<-data_example1
data_example2<-data_example2
data_example3<-data_example3
data_example4<-data_example4
data_example5<-data_example5
data_example6<-data_example6
data_example7<-data_example7
data_example8<-data_example8
data_example9<-data_example9
data_example10<-data_example10
data_example11<-data_example11
data_example12<-data_example12
data_example13<-data_example13
data_example14<-data_example14
data_example15<-data_example15
data_example16<-data_example16
data_example17<-data_example17
data_example18<-data_example18

# Define UI for random distribution app ----
ui <- fluidPage(

  # App title ----
  titlePanel(" "),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # textInput(inputId = " ",
      #            label = " ",
      #          value = " "),

      # Input: Selector for choosing dataset ----
      selectInput("example","Choose example:",
                  choices = c("data_example1", "data_example2", "data_example3",
                              "data_example4","data_example5","data_example6",
                              "data_example7","data_example8","data_example9",
                              "data_example10","data_example11","data_example12",
                              "data_example13","data_example14","data_example15",
                              "data_example16","data_example17","data_example18"


                  )),

      br(),
      varSelectInput("variables", "Choose variable to look at:",example,multiple = TRUE),

      br(),
      varSelectInput("var1","Choose variable in regression:",example),
      br(),
      varSelectInput("var2","Choose variable in regression:",example),
      br(),
      varSelectInput("var3","Choose variable in regression:",example),
      br(),
      sliderInput("timewindow", "Time window:",
                  min = 1, max = 48, value = 7
      ),
      br(),
      numericInput("followuptime", "follow-up Time:",
                   min = 1,max = 3650,value = 84
      ),
      br(),
      selectInput("group", "Arm:",choices = c("treatment", "control"
      )
      )


    ),
    # Main panel for displaying outputs ----
    mainPanel(

      tabsetPanel(type = "tabs",
                  tabPanel("Data", dataTableOutput("table")),
                  tabPanel("MCF plot", plotOutput("MCF")),
                  tabPanel("Table_MCF", dataTableOutput("MCF_table"),dataTableOutput("MCF_diff_table2")),
                  tabPanel("MCF by group", dataTableOutput("MCF_diff_table"),plotOutput("MCF3")),
                  tabPanel("MCF CI estimation", plotOutput("MCF2")),
                  tabPanel("estimation at specific time", dataTableOutput("MCF_specific_time"),
                           dataTableOutput("MCF_specific_time2")))
    )
  )

)




# Define server logic for random distribution app ----


server <- function(input, output,session){


  # Reactive expression to generate the requested distribution ----
  # This is called whenever the inputs change. The output functions
  # defined below then use the value computed from this expression

  d <- reactive({
    switch(input$example,
           "data_example1"=data_example1,
           "data_example2"=data_example2,
           "data_example3"=data_example3,
           "data_example4"=data_example4,
           "data_example5"=data_example5,
           "data_example6"=data_example6,
           "data_example7"=data_example7,
           "data_example8"=data_example8,
           "data_example9"=data_example9,
           "data_example10"=data_example10,
           "data_example11"=data_example11,
           "data_example12"=data_example12,
           "data_example13"=data_example13,
           "data_example14"=data_example14,
           "data_example15"=data_example15,
           "data_example16"=data_example16,
           "data_example17"=data_example17,
           "data_example18"=data_example18


    )

  })

  #dbGetQuery(con, "example",input$example)
  observe ({

    updateSelectInput(session, "variables",
                      choices = colnames(d())
    )

  })

  observe ({

    updateSelectInput(session, "var1",
                      choices = colnames(d())
    )

  })

  observe ({

    updateSelectInput(session, "var2",
                      choices = colnames(d())
    )

  })

  observe ({

    updateSelectInput(session, "var3",
                      choices = colnames(d())
    )

  })






  #output$table <- renderDataTable(d())

  output$table<- renderDataTable({
    if (length(input$variables) == 0) return(d())
    d() %>% dplyr::select(!!!input$variables)
  })


  output$description<-renderDataTable(
    # cbind(c("time1","time2","id" ,"event","terminal","origin"),prepare_descriptive_table(d())$df))
    cbind(names(as.data.frame(d())),prepare_descriptive_table(d())$df))

  output$MCF_table<-renderDataTable({
    item<-d()
    if (is.data.frame(item)){
      item1<-Recur(time=item$days,id=item$id,event=item$event,terminal=item$term_events)~ item$group


    }
    valveMcf0<-mcf(item1)
    valveMcf0@MCF
  })



  output$MCF<-renderPlot({
    item<-d()
    if (is.data.frame(item)){
      item1<-Recur(time=item$days,id=item$id,event=item$event,terminal=item$term_events)~ item$group
      #item$terminal<-rep(0,dim(item)[1])

    }

    # p1<-plot(item1,mcf=T)
    valveMcf0<-mcf(item1)
    df<-as.data.frame(valveMcf0@MCF)
    p1<-ggplot(df, aes(time, MCF)) +
      geom_line(aes(linetype =item.group))+ xlab("Time")+
      theme(axis.text = element_text(size = 10))+
      ggplot2::xlab("Time") + ggplot2::theme_bw()+ theme(axis.text = element_text(size = 20))+
      labs(title="Cummulative incidence rate", y = "event")+
      geom_line(aes(x = time, y =upper,linetype =item.group))+
      geom_line(aes(x = time, y =lower,linetype =item.group))+
      ggtitle("Control vs. Treatment Cummulative incidence rate")+
      theme(axis.text = element_text(size = 10))


    p2<-ggplot(df, aes(time, instRate)) +
      geom_line(aes(linetype =item.group,color=item.group))+ xlab("Time")+
      theme(axis.text = element_text(size = 10))+
      ggplot2::xlab("Time") + ggplot2::theme_bw()+ theme(axis.text = element_text(size = 20))+
      labs(title="Incidence rate", y = "event")+
      #geom_line(aes(x = time, y =upper,linetype =item.group))+
      #geom_line(aes(x = time, y =lower,linetype =item.group))+
      ggtitle("Control vs. Treatment Incidence rate")+
      theme(axis.text = element_text(size = 10))

    df<-df[order(df$item.group,df$time),]
    df1<-df[df$item.group==0,]
    df1<-df1[order(df1$time),]
    df1$instrate.new<-df1$instRate

    timewindow<-input$timewindow
    timepoint.vec<-seq(min(df1$time),max(df1$time),by=timewindow)
    for (j in 1:(length(timepoint.vec)-1)){
      for (i in 1:length(df1$instrate.new)){

        if (df1$time[i]>=timepoint.vec[j] && df1$time[i]<=timepoint.vec[j+1]){
          temp<-df1$instRate[max(which(df1$time<=timepoint.vec[j+1]))]
          df1$instrate.new[i]<-temp

        }

      }
    }

    df2<-df[df$item.group==1,]
    df2<-df2[order(df2$time),]
    df2$instrate.new<-df2$instRate
    #sumup<-function(timepoint,timewindow){
    #  temp<- df$instRate[timepoint:(timepoint+timewindow)]
    #  return(as.numeric(sum(df$instRate[timepoint:(timepoint+timewindow)])))
    #}
    timewindow<-input$timewindow
    #time.cut.point<-round(max(df2$time)/timewindow)
    #sequen<-seq(min(df2$time),max(df2$time),by=time.cut.point)
    timepoint.vec<-seq(min(df2$time),max(df2$time),by=timewindow)
    for (j in 1:(length(timepoint.vec)-1)){
      for (i in 1:length(df2$instrate.new)){

        if (df2$time[i]>=timepoint.vec[j] && df2$time[i]<=timepoint.vec[j+1]){

          temp<-df2$instRate[max(which(df2$time<=timepoint.vec[j+1]))]

          df2$instrate.new[i]<-temp


        }

      }
    }

    df.1<-rbind(df1,df2)
    df.1$instrate.new<-ifelse(df.1$instrate.new<0,0,df.1$instrate.new)


    p3<- ggplot(df.1, aes(time, instrate.new)) +
      geom_line(aes(linetype =item.group,color=item.group))+ xlab("Time")+
      theme(axis.text = element_text(size = 10))+
      ggplot2::xlab("Time") + ggplot2::theme_bw()+ theme(axis.text = element_text(size = 20))+
      labs(title="Incidence rate", y = "event")+
      #geom_line(aes(x = time, y =upper,linetype =item.group))+
      #geom_line(aes(x = time, y =lower,linetype =item.group))+
      ggtitle("Control vs. Treatment Incidence rate")+
      theme(axis.text = element_text(size = 10))


    #p1<-plot(valveMcf0, conf.int = TRUE, col = c("red","royalblue"), lty = c(1, 5)) +
    #  ggtitle("Control vs. Treatment") + xlab("Time")+ theme(axis.text = element_text(size = 10))

    valveMcf<-mcf(item1)
    df<-as.data.frame(valveMcf@MCF)
    p4<-ggplot(df, aes(time, MCF)) +
      geom_line(aes(linetype =item.group))+ xlab("Time")+
      theme(axis.text = element_text(size = 10))+
      ggplot2::xlab("Time") + ggplot2::theme_bw()+ theme(axis.text = element_text(size = 20))+
      labs(title="Cummulative incidence rate", y = "cost")+
      geom_line(aes(x = time, y =upper,linetype =item.group))+
      geom_line(aes(x = time, y =lower,linetype =item.group))+
      ggtitle("Control vs. Treatment Cummulative incidence rate")+
      theme(axis.text = element_text(size = 10))


    p5<-ggplot(df, aes(time, instRate)) +
      geom_line(aes(linetype =item.group,color=item.group))+ xlab("Time")+
      theme(axis.text = element_text(size = 10))+
      ggplot2::xlab("Time") + ggplot2::theme_bw()+ theme(axis.text = element_text(size = 20))+
      labs(title="Incidence rate", y = "cost")+
      #geom_line(aes(x = time, y =upper,linetype =item.group))+
      #geom_line(aes(x = time, y =lower,linetype =item.group))+
      ggtitle("Control vs. Treatment Incidence rate")+
      theme(axis.text = element_text(size = 10))



    df<-df[order(df$item.group,df$time),]
    df1<-df[df$item.group==0,]
    df1<-df1[order(df1$time),]
    df1$instrate.new<-df1$instRate
    #df1$instrate.new<-rep(0,dim(df1)[1])
    #sumup<-function(timepoint,timewindow){
    #  temp<- df$instRate[timepoint:(timepoint+timewindow)]
    #  return(as.numeric(sum(df$instRate[timepoint:(timepoint+timewindow)])))
    #}
    timewindow<-input$timewindow
    timepoint.vec<-seq(min(df1$time),max(df1$time),by=timewindow)
    for (j in 1:(length(timepoint.vec)-1)){
      for (i in 1:length(df1$instrate.new)){


        if (df1$time[i]>=timepoint.vec[j] && df1$time[i]<=timepoint.vec[j+1]){

          temp<-df1$instRate[max(which(df1$time<=timepoint.vec[j+1]))]
          df1$instrate.new[i]<-temp

        }

      }

    }

    df2<-df[df$item.group==1,]
    df2<-df2[order(df2$time),]

    df2$instrate.new<-df2$instRate
    #df2$instrate.new<-rep(0,dim(df2)[1])
    #sumup<-function(timepoint,timewindow){
    #  temp<- df$instRate[timepoint:(timepoint+timewindow)]
    #  return(as.numeric(sum(df$instRate[timepoint:(timepoint+timewindow)])))
    #}
    timewindow<-input$timewindow
    timepoint.vec<-seq(min(df2$time),max(df2$time),by=timewindow)
    for (j in 1:(length(timepoint.vec)-1))  {
      for (i in 1:length(df2$instrate.new)){

        if (df2$time[i]>=timepoint.vec[j] &&  df2$time[i]<=timepoint.vec[j+1]){

          temp<-df2$instRate[max(which(df2$time<=timepoint.vec[j+1]))]
          df2$instrate.new[i]<-temp

        }


      }

    }


    df.2<-rbind(df1,df2)
    df.2$instrate.new<-ifelse(df.2$instrate.new<0,0,df.2$instrate.new)


    p6<- ggplot(df.2, aes(time, instrate.new)) +
      geom_line(aes(linetype =item.group,color=item.group))+ xlab("Time")+
      theme(axis.text = element_text(size = 10))+
      ggplot2::xlab("Time") + ggplot2::theme_bw()+ theme(axis.text = element_text(size = 20))+
      labs(title="Incidence rate", y = "cost")+
      #geom_line(aes(x = time, y =upper,linetype =item.group))+
      #geom_line(aes(x = time, y =lower,linetype =item.group))+
      ggtitle("Control vs. Treatment Incidence rate")+
      theme(axis.text = element_text(size = 10))

    df.1$ratio.se<-0
    df.1$ratio<-0
    df.1<-df.1[order(df.1$item.group,df.1$time),]
    vec1<-df.1[df.1$item.group==0,]$instrate.new
    vec2<-df.1[df.1$item.group==1,]$instrate.new
    vec3<-df.1[df.1$item.group==0,]$se
    vec4<-df.1[df.1$item.group==1,]$se
    tempt<-vector()
    tempt2<-vector()
    for (i in 1:min(length(vec1),length(vec2))){


      tempt[i]<-ifelse(vec1[i]!=0,vec2[i]/vec1[i],0)
      tempt2[i]<-ifelse(vec1[i]!=0,vec2[i]/vec1[i]*sqrt(vec3[i]^2/vec1[i]^2+vec4[i]^2/vec2[i]^2),0)
    }

    df.1$ratio[1:length(tempt)]<-tempt
    df.1$ratio.se[1:length(tempt2)]<-tempt2
    df.1$ratio.up<-df.1$ratio+1.96*df.1$ratio.se
    df.1$ratio.lw<-df.1$ratio-1.96*df.1$ratio.se


    df.2$ratio.se<-0
    df.2$ratio<-0
    df.2<-df.2[order(df.2$item.group,df.2$time),]
    vec1<-df.2[df.2$item.group==0,]$instrate.new
    vec2<-df.2[df.2$item.group==1,]$instrate.new
    vec3<-df.2[df.2$item.group==0,]$se
    vec4<-df.2[df.2$item.group==1,]$se
    tempt<-vector()
    tempt2<-vector()
    for (i in 1:min(length(vec1),length(vec2))){


      tempt[i]<-ifelse(vec1[i]!=0,vec2[i]/vec1[i],0)
      tempt2[i]<-ifelse(vec1[i]!=0,vec2[i]/vec1[i]*sqrt(vec3[i]^2/vec1[i]^2+vec4[i]^2/vec2[i]^2),0)

    }
    df.2$ratio[1:length(tempt)]<-tempt
    df.2$ratio.se[1:length(tempt2)]<-tempt2
    df.2$ratio.up<-df.2$ratio+1.96*df.2$ratio.se
    df.2$ratio.lw<-df.2$ratio-1.96*df.2$ratio.se


    p7<- ggplot(df.1, aes(time, ratio)) +
      geom_line(color="red")+ xlab("Time")+
      #geom_line(aes(time, ratio.up,linetype ="dashed"))+
      #geom_line(aes(time, ratio.lw,linetype ="dashed"))+
      theme(axis.text = element_text(size = 10))+
      ggplot2::xlab("Time") + ggplot2::theme_bw()+ theme(axis.text = element_text(size = 20))+
      labs(title="Incidence rate ratio", y = "Visit")+
      #geom_line(aes(x = time, y =upper,linetype =item.group))+
      #geom_line(aes(x = time, y =lower,linetype =item.group))+
      ggtitle("Control vs. Treatment Incidence rate ratio")+
      theme(axis.text = element_text(size = 10))

    p8<- ggplot(df.2, aes(time, ratio)) +
      geom_line(color="red")+ xlab("Time")+
      #geom_line(aes(time, ratio.up,linetype ="dashed"))+
      #geom_line(aes(time, ratio.lw,linetype ="dashed"))+
      theme(axis.text = element_text(size = 10))+
      ggplot2::xlab("Time") + ggplot2::theme_bw()+ theme(axis.text = element_text(size = 20))+
      labs(title="Incidence rate ratio", y = "cost")+
      #geom_line(aes(x = time, y =upper,linetype =item.group))+
      #geom_line(aes(x = time, y =lower,linetype =item.group))+
      ggtitle("Control vs. Treatment Incidence rate ratio")+
      theme(axis.text = element_text(size = 10))

    ggarrange(p1,p2,p3,p7,ncol = 2, nrow = 2)

  }, height = 800, width = 800 )




  output$Ratio_visit<-renderDataTable({
    item<-d()
    if (is.data.frame(item)){
      item1<-Recur(time=item$days,id=item$id,event=item$event,terminal=item$term_events)~ item$group

    }

    # p1<-plot(item1,mcf=T)
    valveMcf0<-mcf(item1)
    df<-as.data.frame(valveMcf0@MCF)

    df<-df[order(df$item.group,df$time),]
    df1<-df[df$item.group==0,]
    df1<-df1[order(df1$time),]
    df1$instrate.new<-df1$instRate

    timewindow<-input$timewindow
    timepoint.vec<-seq(min(df1$time),max(df1$time),by=timewindow)
    for (j in 1:(length(timepoint.vec)-1)){
      for (i in 1:length(df1$instrate.new)){

        if (df1$time[i]>=timepoint.vec[j] && df1$time[i]<=timepoint.vec[j+1]){
          temp<-df1$instRate[max(which(df1$time<=timepoint.vec[j+1]))]
          df1$instrate.new[i]<-temp

        }

      }
    }

    df2<-df[df$item.group==1,]
    df2<-df2[order(df2$time),]
    df2$instrate.new<-df2$instRate
    #sumup<-function(timepoint,timewindow){
    #  temp<- df$instRate[timepoint:(timepoint+timewindow)]
    #  return(as.numeric(sum(df$instRate[timepoint:(timepoint+timewindow)])))
    #}
    timewindow<-input$timewindow
    #time.cut.point<-round(max(df2$time)/timewindow)
    #sequen<-seq(min(df2$time),max(df2$time),by=time.cut.point)
    timepoint.vec<-seq(min(df2$time),max(df2$time),by=timewindow)
    for (j in 1:(length(timepoint.vec)-1)){
      for (i in 1:length(df2$instrate.new)){

        if (df2$time[i]>=timepoint.vec[j] && df2$time[i]<=timepoint.vec[j+1]){

          temp<-df2$instRate[max(which(df2$time<=timepoint.vec[j+1]))]

          df2$instrate.new[i]<-temp


        }

      }
    }

    df.1<-rbind(df1,df2)
    df.1$instrate.new<-ifelse(df.1$instrate.new<0,0,df.1$instrate.new)


    df.1$ratio.se<-0
    df.1$ratio<-0
    df.1<-df.1[order(df.1$item.group,df.1$time),]
    vec1<-df.1[df.1$item.group==0,]$instrate.new
    vec2<-df.1[df.1$item.group==1,]$instrate.new
    vec3<-df.1[df.1$item.group==0,]$se
    vec4<-df.1[df.1$item.group==1,]$se
    tempt<-vector()
    tempt2<-vector()
    for (i in 1:min(length(vec1),length(vec2))){


      tempt[i]<-ifelse(vec1[i]!=0,vec2[i]/vec1[i],0)
      tempt2[i]<-ifelse(vec1[i]!=0,vec2[i]/vec1[i]*sqrt(vec3[i]^2/vec1[i]^2+vec4[i]^2/vec2[i]^2),0)
    }

    df.1$ratio[1:length(tempt)]<-tempt
    df.1$ratio.se[1:length(tempt2)]<-tempt2
    df.1$ratio.up<-df.1$ratio+1.96*df.1$ratio.se
    df.1$ratio.lw<-df.1$ratio-1.96*df.1$ratio.se
    res<-as.data.frame(df.1$ratio)
    colnames(res)<-"Visit.Ratio"
    res
  })




  output$MCF2<-renderPlot({
    item<-d()
    if (is.data.frame(item)){
      item<-Recur(time=item$days,id=item$id,event=item$event,terminal=item$term_events)
    }
    temp<-mcf(item~1)
    temp2<-mcf(item ~ 1,variance = "Poisson")
    temp3<-mcf(item ~ 1, variance = "bootstrap", control = list(B = 1e3))
    ## comparing the standard error estimates

    ciDat <- rbind(cbind(temp@MCF, Method = "Lawless & Nadeau"),
                   cbind(temp2@MCF, Method = "Poisson"),
                   cbind(temp3@MCF, Method = "Bootstrap"))
    p1<- ggplot(ciDat, aes(x = time, y = se)) +
      geom_step(aes(color = Method, linetype = Method)) +
      xlab("Time") + ylab("SE estimates") + theme_bw()+ theme(axis.text = element_text(size = 20))

    p2<-ggplot(ciDat, aes(x = time)) +
      geom_step(aes(y = MCF), color = "grey") +
      geom_step(aes(y = lower, color = Method, linetype = Method)) +
      geom_step(aes(y = upper, color = Method, linetype = Method)) +
      xlab("Time") + ylab("Confidence intervals") + theme_bw()+ theme(axis.text = element_text(size = 20))
    ggarrange(p1,ncol = 1, nrow = 1)
  }, height = 400, width = 600 )


  output$MCF3<-renderPlot({
    b02_tr<-d()
    simuMcf <- mcf(Recur(days,id, event) ~ group ,
                   data = b02_tr)
    p1<-plot(simuMcf, conf.int = TRUE, lty = 1:4, legendName = "Treatment")+ theme(axis.text = element_text(size = 10))
    ## one sample MCF object of two groups
    mcf0 <- mcf(Recur(days, id, event) ~ group, data = b02_tr)
    mcf_diff0 <- mcfDiff(mcf0)
    p2<-plot(mcf_diff0)+ theme(axis.text = element_text(size = 10))




    #constFit <- rateReg(Recur(days, id, event) ~ as.factor(group) , data = b02_tr)
    #twoPiecesFit <- rateReg(Recur(time1, id, event) ~ group, df = 2,
    #                       data = b02_tr)

    #piecesFit <- rateReg(Recur(days, id, event) ~ group, data = b02_tr)
    #splineFit <- rateReg(Recur(days, id, event) ~ group , data = b02_tr)

    #newDat <- data.frame(x1 = c(0, 0), group = c("Treat","Contr"))

    #estmcf <- mcf(splineFit, newdata = newDat, groupName = "Group",
    #              groupLevels = c("Treatment", "Control"))
    #p4<-plot(estmcf, conf.int = TRUE, col = c("royalblue", "red"), lty = c(1, 5)) +
    #  ggtitle("Control vs. Treatment using B-spine") + xlab("Time")+ theme(axis.text = element_text(size = 10))

    #p5<-plot(estmcf) +
    #  geom_ribbon(data = estmcf@MCF, alpha = 0.2,
    #              aes(x = time, ymin = lower, ymax = upper, fill = Group)) +
    #  ggtitle("Control vs. Treatment using B-spine") + xlab("Time")+ theme(axis.text = element_text(size = 10))

    #df<-estmcf@MCF
    #df$insrate<-c(rep(0,dim(df)[1]))
    #df<-df[order(df$Group,df$time),]
    #for (i in 2:length(df$insrate)){
    #  df$insrate[i]<-df$MCF[i]-df$MCF[i-1]
    #}
    #df$insrate<-ifelse(df$insrate<0,0,df$insrate)



    #p6<-ggplot(df, aes(time, insrate)) +
    #  geom_line(aes(linetype =Group,color=Group))+ xlab("Time")+
    # theme(axis.text = element_text(size = 10))+
    #  ggplot2::xlab("Time") +
    #  labs(title="Incidence rate", y = "Incidence rate")+
    #geom_line(aes(x = time, y =upper,linetype =Group))+
    #geom_line(aes(x = time, y =lower,linetype =Group))+
    #    ggtitle("Control vs. Treatment")+
    #    theme(axis.text = element_text(size = 10))

    ggarrange(p1,p2,
              ncol = 1, nrow = 2)
  }, height = 800, width = 800 )
  # Generate a plot of the data ----
  # Also uses the inputs to build the plot label. Note that the
  # dependencies on the inputs and the data reactive expression are
  # both tracked, and all expressions are called in the sequence
  # implied by the dependency graph.

  output$MCF_diff_table<-renderDataTable({
    b02_tr<-d()
    mcf0 <- mcf(Recur(days, id, event) ~ group, data = b02_tr)
    mcf_diff0 <- mcfDiff(mcf0)
    round(mcf_diff0@test,4)
  })

  output$MCF_diff_table2<-renderDataTable({
    b02_tr<-d()
    mcf0 <- mcf(Recur(days, id, event) ~ group, data = b02_tr)
    mcf_diff0 <- mcfDiff(mcf0)
    mcf_diff0@MCF
  })

  # output$effect<-renderDataTable({
  #  b02_tr<-d()
  #splineFit <- rateReg(Recur(time2, id, event) ~ group , data = b02_tr,
  #                    spline = "bSp", degree = 3L, knots = seq(from = 28, to = 140, by = 28))
  #  splineFit <- rateReg(Recur(days, id, event) ~ group , data = b02_tr)

  #b02_tr$cost<-ifelse(b02_tr$event==1 & b02_tr$group==0,4,b02_tr$cost)

  #  round(splineFit@estimates$beta,4)
  #})

  output$MCF_specific_time<-renderDataTable({
    b02_tr<-d()
    mcf0 <- mcf(Recur(days,id,event) ~ as.factor(group), data = b02_tr)
    mcf_diff0 <- mcfDiff(mcf0)
    df<-mcf_diff0@MCF
    colnames(df)<-c("time","difference in MCF","se","lower","upper")
    time_diff<-as.numeric(input$followuptime-df$time)
    time<-df$time[which(time_diff>=0)]
    #df[df$time==head(time)[1],]
    #round(df[df$time==head(time)[1],],4)
    round(df[df$time==max(time),],4)
  })


  output$MCF_specific_time2<-renderDataTable({
    b02_tr<-d()
    mcf0 <- mcf(Recur(days,id,event) ~ as.factor(group), data = b02_tr)
    df<-mcf0@MCF
    colnames(df)<-c("time","numRisk","instRate","MCF","se","lower","upper", "as.factor.group.")
    df<-df[df$as.factor.group.==ifelse(input$group=="treatment",1,0),]
    time_diff<-as.numeric(input$followuptime-df$time)
    time<-df$time[which(time_diff>=0)]
    df$as.factor.group.<-as.numeric(df$as.factor.group.)
    #df[df$time==head(time)[1],]
    #round(df[df$time==head(time)[1],],4)
    round(df[df$time==max(time),],4)
  })



}



# Run the application
shinyApp(ui = ui, server = server)

