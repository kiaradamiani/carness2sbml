__author__ = 'AirZoi'
# -*- coding: cp1252 -*-
import os
import csv
import numpy
import tarfile



def getLastReactionID(lastGen,filesPath):

    #apri file times ultima generazione per sapere l'id dell'ultima reazione
    file_times="%s%s%s%s" % (filesPath,'times_',lastGen,'_1.csv')
    reactions_id,time=numpy.loadtxt(file_times, dtype=int, delimiter='\t', usecols=(0,1), unpack=True)
    id_last_reaction=reactions_id[len(reactions_id)-1]+1

    print 'id_last_reaction',id_last_reaction

    #conta cifre dell'id reazione per sapere quanti zeri ci sono nel file
    number=id_last_reaction
    n_cifre=0
    while number>=1:
        number=number/10
        n_cifre=n_cifre+1
        #aggiungo zero per le cifre che rimangono fino a 9
    zeros=''
    for e1 in range(0,9-n_cifre):
        zeros="%s%s" % (zeros,'0')

    return "%s%s" %(zeros,id_last_reaction)


def getParam(paramName,file_conf):

    #leggo file conf
    # ##########################################
    #
    for line in open(file_conf,'r').readlines():
        thisline = line.split("=");

        if thisline[0]==paramName:
            return thisline[1].strip()


def importChemistry(lastGen,filesPath=""):

    id_last_reaction=getLastReactionID(lastGen,filesPath)


    file_conf=filesPath+'acsm2s.conf'
    lunghezza_max_fd=getParam('lunghezza_max_fd',file_conf)
    InitialNumberOfSpecies=pow(2,int(lunghezza_max_fd)+1)-2

    overallConcentration=getParam('overallConcentration',file_conf)
    initialConc=float(overallConcentration)/float(InitialNumberOfSpecies)

    #nome file spesis ultima generazione e ultima reazione
    file_species="%s%s%s%s%s%s" % (filesPath,'species_',lastGen,'_1_',id_last_reaction,'.csv')
    file_reactions="%s%s%s%s%s%s" % (filesPath,'reactions_',lastGen,'_1_',id_last_reaction,'.csv')
    file_catalysis="%s%s%s%s%s%s" % (filesPath,'catalysis_',lastGen,'_1_',id_last_reaction,'.csv')
    file_conf="%s%s" %(filesPath,'acsm2s.conf')
    file_influx="%s%s" %(filesPath,'_acsinflux.conf')
    file2write='ExportFile.dizzy'

    with open(file2write,'w') as file2write:

        comment= "\n"+ "\n"+'//SPECIES'+ "\n"
        file2write.write(comment)

        print comment

        #costruisco dizionario con id specie e nome specie (e dizionario con nome specie e indicazione se bufferizzata) e scrivo elenco specie con concentrazioni su file
        speciesIDs={}
        speciesNames={}

        for line in open(file_species,'r').readlines():
            thisline = line.split("\t");

            speciesID=thisline[0]
            speciesName=thisline[1]
            lockedConcentration=thisline[14]
            Kdecomp=thisline[5]


#            #se la concentrazione è bloccata la specie deve comparire nelle reazioni con il dollaro davanti
#            if (int(lockedConcentration) == 1):
#                speciesIDs[speciesID]='$'+speciesName
#            else:
#                speciesIDs[speciesID]=speciesName

            speciesIDs[speciesID]=speciesName
            speciesNames[speciesName]=lockedConcentration.strip()


            #se specie non era presente all'inzio metto concetrazione a 0
            if int(speciesID)<int(InitialNumberOfSpecies):
                initialConcentration=thisline[1]+'='+str(initialConc)+';'
            else:
                initialConcentration=thisline[1]+'='+str(0)+';'

            print initialConcentration
            file2write.write(initialConcentration)
            file2write.write('\n')


        #costruisco dizionario con id reaction e specie catalyst(già tradotta dall'id) e rate
        catalysisIDs={}
        for line in open(file_catalysis,'r').readlines():
            thisline = line.split("\t");
            catalysisIDs[thisline[2]]= speciesIDs[thisline[1]]
            Kcleav=thisline[5]
            Kcond=thisline[4]
            Kcomp=thisline[6].strip()


        comment=  "\n"+ "\n"+'//reactions'+ "\n"
        print comment
        file2write.write(comment)


        #leggo reazioni
        line_counter=0
        for line in open(file_reactions,'r').readlines():
            thisline = line.split("\t");
            reactionType=thisline[1]

            #se è un cleavage
            if reactionType=='1':
                products="" #inizializzo prodotti

                reactionID=thisline[0]
                substrate=speciesIDs[thisline[2]]
                #se concentrazione è bloccata aggiungo dollaro
                if speciesNames[substrate]=='1':
                    substate='$'+substrate

                product1=speciesIDs[thisline[3]]
                print "//Cleavage product 1 ",product1
                if (speciesNames[product1]=='1'):
                    product1='$'+product1
                else:
                    products=product1

                product2=speciesIDs[thisline[4]]
                print "//Cleavage product 2 ",product2

                if speciesNames[product2]=='1':
                    product2='$'+product2
                else:
                    if len(products)>0:
                        products=products+"+"+product2
                    else:
                        products=product2

                catalyst=catalysisIDs[reactionID]
                if (speciesNames[catalyst]=='1'):
                    catalyst='$'+catalyst
                else:
                    if len(products)>0:
                        products=products+"+"+catalyst
                    else:
                        products=catalyst

                cleavage='cleavage'+reactionID+','+substrate+'+'+catalyst+' -> '+products+', '+Kcleav+';'

                file2write.write(cleavage)
                file2write.write('\n')

                print cleavage

            #se è una condensation
            products=""
            decomplex=""
            if reactionType=='0':
                reactionID=thisline[0]
                substrate1=speciesIDs[thisline[3]]
                if speciesNames[substrate1]=='1':
                    substrate1='$'+substrate1
                else:
                    decomplex=substrate1

                substrate2=speciesIDs[thisline[4]]
                if speciesNames[substrate2]=='1':
                    substrate2='$'+substrate2

                product=speciesIDs[thisline[2]]

                print '//product condendation',product

                if product not in speciesNames.keys():
                    speciesNames[product]='0'
                    productConc=product+"="+0+";"
                    print productConc
                    file2write.write(productConc)
                    products=product
                else:
                    if speciesNames[product]=='1':
                        print 'speciesNames[product]==1:'
                        product='$'+product
                    else:
                        products=product

                catalyst=catalysisIDs[reactionID]
                if speciesNames[catalyst]=='1':
                    catalyst='$'+catalyst
                else:
                    if len(products)>0:
                        products=product+"+"+catalyst
                    else:
                        products=catalyst
                    if len(decomplex)>0:
                        decomplex=decomplex+"+"+catalyst
                    else:
                        decomplex=catalyst


                complex=catalyst+substrate1.strip('$')
                #se complesso non è già nel file species devo settare concentrazione iniziale
                if complex not in speciesNames.keys():
                    speciesNames[complex]='0'
                    complexConc=complex+"="+'0'+";"+"\n"
                    print complexConc
                    file2write.write(complexConc)

                complexation= 'complexation'+reactionID+','+substrate1+'+'+catalyst+' ->'+complex+','+Kcomp+';'
                decomplexation='decomplexation'+reactionID+','+complex+' ->'+decomplex+','+Kdecomp+';'
                condensation='condensation'+reactionID+','+complex+'+'+substrate2+  ' ->' +products+','+Kcond+';'

                print complexation
                file2write.write(complexation)
                file2write.write('\n')

                print decomplexation
                file2write.write(decomplexation)
                file2write.write('\n')

                print condensation
                file2write.write(condensation)
                file2write.write('\n')




                #with open(file2write,'w') as file2write:

                #file2write.write(substrate1,'+',catalyst,' ->', catalyst.strip(),substrate1.strip(),',',Kcomp,';')

            #print catalyst,substrate1,' ->',substrate1,'+',catalyst,',',Kdecomp,';'
            #print catalyst,substrate1,'+',substrate2,  ' ->' ,substrate1,substrate2,'+' ,catalyst,',',Kcond,';'


#    with open(file2write,'w') as file2write:
#
#    #leggo file REACTIONS
#    ##    #########################################
#    #creo vettore delle reazioni in cui mettere se sono cleavage o cond
#
#        reactions_id_type,reactions_type=numpy.loadtxt(file_reactions, dtype=int, delimiter='\t', usecols=(0,1), unpack=True)
#
#        line_counter=0
#        for line in open(file_reactions,'r').readlines():
#            line_counter=line_counter+1
#            thisline = line.split("\t");
#
#            thisline[5]=0 #azzero contatore
#            #scrivo riga su file
#            for col in range(0,len(thisline)):
#                file2write.write(str(thisline[col]))
#                if col<len(thisline)-1: #se non Ë ultima colonna metto tab
#                    file2write.write('\t')



importChemistry(100)

#
#    ##########################################
#
#    #leggo file species vecchi e scrivo il nuovo
#    #########################################
#    percorso_file_species_nuovo="%s%s" % (nome_nuova_rete,'\\_acsspecies.csv')
#    file2write = open(percorso_file_species_nuovo, "w")
#
#    line_counter=0
#    for line in open(percorso_file_species,'r').readlines():
#        line_counter=line_counter+1
#        thisline = line.split("\t");
#
#        #per le prime 29 righe cambio valore colonna concentrazione (2) a .333333e-0
#        if line_counter<31:
#            thisline[2]=1.110000e-03 #cambio valore colonna concentrazione (2) a .333333e-0
#        #per le altre righe lo cambio a 0
#        else:
#            thisline[2]=0
#
#        thisline[8]=0 #metto species age [colonna 8] =0
#        thisline[5]=K_cpx #metto K_cpx uguale al parametro che voglio
#
#        #scrivo riga su file
#        for col in range(0,len(thisline)):
#            file2write.write(str(thisline[col]))
#            if col<len(thisline)-1: #se non Ë ultima colonna metto tab
#                file2write.write('\t')
#    file2write.close()
#    ##########################################
#
#
#    #leggo file catalysis vecchi e scrivo il nuovo
#    #########################################
#
#    #nome file spesis ultima generazione e ultima reazione
#    nome_file_catalysis="%s%s%s%s%s%s" % ('\\res\\catalysis_',lastGen,'_1_',zeros,id_last_reaction,'.csv')
#    percorso_file_catalysis="%s%s" % (nome_rete, nome_file_catalysis)
#    print percorso_file_catalysis
#
#    percorso_file_catalysis_nuovo="%s%s" % (nome_nuova_rete,'\\_acscatalysis.csv')
#    file2write = open(percorso_file_catalysis_nuovo, "w")
#
#    line_counter=0
#    for line in open(percorso_file_catalysis,'r').readlines():
#        line_counter=line_counter+1
#        thisline = line.split("\t");
#
#        idr=int(thisline[2]) #id reazione catalizzata
#        idr_type=reactions_type[idr]
#
#        #se la reazione Ë una condensazione
#        if (idr_type==0):
#            thisline[5]=float(Kdiss)/float(revRctRatio)
#
#        if (idr_type==1):
#            thisline[4]=float(Kass)/float(revRctRatio)
#            thisline[6]="%d%s" % (float(Kcpx)/float(revRctRatio), '\n') #risetto Kcpx (siccome Ë ultima colonna devo anche mettere a capo)
#
#        thisline[3]=0 #azzero contatore catalisi
#
#        #scrivo riga su file
#        for col in range(0,len(thisline)):
#            file2write.write(str(thisline[col]))
#            if col<len(thisline)-1: #se non Ë ultima colonna metto tab
#                file2write.write('\t')
#    file2write.close()
#    ##########################################
#
#
#
#    #modifico file CONF
#    ##########################################
#
#    nome_file_conf='\\acsm2s.conf'
#    percorso_file_conf="%s%s" % (nome_rete, nome_file_conf)
#    percorso_file_conf_nuovo="%s%s" % (nome_nuova_rete,'\\acsm2s.conf')
#
#    file2write = open(percorso_file_conf_nuovo, "w")
#
#    line_counter=0
#    for line in open(percorso_file_conf,'r').readlines():
#        line_counter=line_counter+1
#        thisline = line.split("=");
#
#        if thisline[0]=='K_ass':
#            thisline[1]="%s%s" % (Kass, '\n')
#        if thisline[0]=='K_diss':
#            thisline[1]="%s%s" % (Kdiss, '\n')
#        if thisline[0]=='K_cpx':
#            thisline[1]="%s%s" % (Kcpx, '\n')
#        if thisline[0]=='K_cpxDiss':
#            thisline[1]="%s%s" % (K_cpx, '\n')
#        if thisline[0]=='revRctRatio':
#            thisline[1]="%s%s" % (revRctRatio, '\n')
#
#
#        #scrivo riga su file
#        for col in range(0,len(thisline)):
#            file2write.write(str(thisline[col]))
#            if col<len(thisline)-1: #se non Ë ultima colonna metto tab
#                file2write.write('=')
#    file2write.close()
#
#    ##########################################
#
#
