__author__ = 'AirZoi'
# -*- coding: utf-8 -*-


import os
import csv
import numpy
import tarfile
from importChemistryParameters import *


#da file times recuper Id ultima reazione
def getLastReactionID(lastGen,ResFilesPath):

    #apri file times ultima generazione per sapere l'id dell'ultima reazione
    file_times="%s%s%s%s" % (ResFilesPath,'times_',lastGen,'_1.csv')
    reactions_id,time=numpy.loadtxt(file_times, dtype=int, delimiter='\t', usecols=(0,1), unpack=True)
    id_last_reaction=reactions_id[len(reactions_id)-1]

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

#Da file conf recupera valore parametro che si vuole
def getParam(paramName,file_conf):

    #leggo file conf
    # ##########################################
    #
    for line in open(file_conf,'r').readlines():
        thisline = line.split("=");

        if thisline[0]==paramName:
            return thisline[1].strip()

#leggo chimica e creo file dizzy
def importChemistry(lastGen,filesPath,ConfFilePath,outputFilePath):
    #id ultima reazione
    id_last_reaction=getLastReactionID(lastGen,filesPath)

    #nome files e parametri
    #######################
    file_conf=ConfFilePath+'acsm2s.conf'
    lunghezza_max_fd=getParam('lunghezza_max_fd',file_conf)
    V1=float(getParam('volume',file_conf))
    InitialNumberOfSpecies=pow(2,int(lunghezza_max_fd)+1)-2
    
    overallConcentration=getParam('overallConcentration',file_conf)
    initialConc=(float(overallConcentration)/float(InitialNumberOfSpecies))*float(V1) #moltiplico concentrazione per volume perchè SBML vuole quantità in AMOUNT

    #nome file spesis ultima generazione e ultima reazione
    file_species="%s%s%s%s%s%s" % (filesPath,'species_',lastGen,'_1_',id_last_reaction,'.csv')
    file_reactions="%s%s%s%s%s%s" % (filesPath,'reactions_',lastGen,'_1_',id_last_reaction,'.csv')
    file_catalysis="%s%s%s%s%s%s" % (filesPath,'catalysis_',lastGen,'_1_',id_last_reaction,'.csv')
    file_conf="%s%s" %(filesPath,'acsm2s.conf')
    file_influx="%s%s" %(filesPath,'_acsinflux.conf')
    
    
    #inizio a scrivere reazioni su file dizzy
    file2write=outputFilePath+'ExportFile.dizzy'

    with open(file2write,'w') as file2write:
        
        #scrivo volume compartimento
        volume= "V1 = "+ str(float(V1))+ ";\n"
        print volume
        file2write.write(volume)

        comment= "\n"+ "\n"+'//SPECIES'+ "\n"
        file2write.write(comment)

        print comment

        #costruisco dizionario con id specie e nome specie (e dizionario con nome specie e indicazione se bufferizzata) e scrivo elenco specie con concentrazioni su file
        speciesIDs={}
        speciesNames={}
        writtenComparments={} #siccome non posso specificare compartimento pi˘ di una volta ne tengo traccia

        #per ogni specie nel file scies
        for line in open(file_species,'r').readlines():
            thisline = line.split("\t");

            speciesID=thisline[0]
            speciesName=thisline[1]
            lockedConcentration=thisline[14]
            Kdecomp=str(float(thisline[5])/V1)

            #aggiung nome specie a dizionario id - nomi
            speciesIDs[speciesID]=speciesName
            #aggiungo se è bufferizzato no al dizionario con i nomi
            speciesNames[speciesName]=lockedConcentration.strip()

            print 'speciesNames[speciesName] ',speciesName


            #se specie fa parte del firing disk la concentrazione iniziale sarà  overallConcentration/numberofSpecies altrimenti sarà 0
            if int(speciesID)<int(InitialNumberOfSpecies):
                initialConcentration=thisline[1]+'='+str(initialConc)+';'
            else:
                initialConcentration=thisline[1]+'='+str(0)+';'


            #se non ho già scritto scrivo quanittà iniziale specie e compartimento di riferimento
            if thisline[1] not in writtenComparments.keys():
                print initialConcentration
                file2write.write(initialConcentration)
                file2write.write('\n')
                
                print " if thisline[1] not in writtenComparments.keys():"

                speciesCompartment=thisline[1]+" @ V1;"+"\n"
                print speciesCompartment
                writtenComparments[thisline[1]]='V1'
                file2write.write(speciesCompartment)
                
                file2write.write('//thisline[1] not in writtenComparments.keys(): \n')


        #costruisco dizionario con id reaction e nome specie catalyst(gi‡ tradotta dall'id) e rate
        catalysisIDs={}
        for line in open(file_catalysis,'r').readlines():
            thisline = line.split("\t");
            catalysisIDs[thisline[2]]= speciesIDs[thisline[1]]
            #recupero costanti cinetiche
            Kcleav=str(float(thisline[5])/V1)
            Kcond=str(float(thisline[4])/V1)
            Kcomp=str(float(thisline[6].strip())/V1)


        comment=  "\n"+ "\n"+'//reactions'+ "\n"
        print comment
        file2write.write(comment)


        #leggo reazioni (per ogni reazione=
        line_counter=0
        for line in open(file_reactions,'r').readlines():
            thisline = line.split("\t");
            reactionType=thisline[1]

            #se Ë un cleavage
            if reactionType=='1':
                products="" #inizializzo prodotti

                reactionID=thisline[0]
                substrate=speciesIDs[thisline[2]]
                #se concentrazione Ë bloccata aggiungo dollaro
                if speciesNames[substrate]=='1':
                    ##print "==1"
                    substrate='$'+substrate

                product1=speciesIDs[thisline[3]]
                print "//Cleavage product 1 ",product1
                if (speciesNames[product1.strip('$')]=='1'):
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
                #print 'catalyst',catalyst
                #print 'speciesNames[catalyst]',speciesNames[catalyst]
                if (speciesNames[catalyst]=='1'):
                    print "==1"
                    #catalyst='$'+catalyst
                else:
                    if len(products)>0:
                        products=products+"+"+catalyst
                    else:
                        products=catalyst

                cleavage='cleavage'+reactionID+','+substrate+'+'+catalyst+' -> '+products+', '+Kcleav+';'

                file2write.write(cleavage)
                file2write.write('\n')

                print cleavage

            #se Ë una condensation
            products=""
            decomplex=""
            if reactionType=='0':
                reactionID=thisline[0]
                substrate1=speciesIDs[thisline[3]]
                if speciesNames[substrate1.strip('$')]=='1':
                    substrate1='$'+substrate1
                else:
                    decomplex=substrate1

                substrate2=speciesIDs[thisline[4]]
                if speciesNames[substrate2.strip('$')]=='1':
                    substrate2='$'+substrate2

                product=speciesIDs[thisline[2]]

                print '//product condendation',product
                
                #se non ho già settato cocentrzione e compartment prodotto lo faccio
                if product not in writtenComparments.keys():
                    print "if product not in speciesNames.keys():"

                    speciesNames[product.strip('$')]='0'
                    productConc=product+"="+0+";"
                    print productConc
                    file2write.write(productConc)
                    products=product
                    speciesCompartment=product+" @ V1;"+"\n"
                    file2write.write(speciesCompartment)
                    print speciesCompartment
                    
                    #aggiungo a dizionario delle specie già scritte
                    writtenComparments[product.strip('$')]='V1'

                else:
                    if speciesNames[product.strip('$')]=='1':
                        product='$'+product
                    else:
                        products=product

                catalyst=catalysisIDs[reactionID]
                if speciesNames[catalyst.strip('$')]=='1':
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


                complex=catalyst.strip('$')+substrate1.strip('$')
                #se complesso non Ë gi‡ nel file species devo settare concentrazione iniziale
                if complex.strip('$') not in writtenComparments.keys():
                    print 'if complex not in speciesNames.keys():'

                    speciesNames[complex.strip('$')]='0'
                    complexConc=complex.strip('$')+"="+'0'+";"+"\n"
                    print complexConc
                    file2write.write(complexConc)
                    speciesCompartment=complex.strip('$')+" @ V1;"+"\n"
                    file2write.write(speciesCompartment)
                    print speciesCompartment
                    
                    #aggiungo a dizionario delle specie già scritte
                    writtenComparments[complex.strip('$')]='V1'



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




importChemistry(lastGen,ResFilesPath,ConfFilePath,outputFilePath)
