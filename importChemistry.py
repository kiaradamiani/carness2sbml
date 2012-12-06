__author__ = 'AirZoi'
# -*- coding: cp1252 -*-
import os
import csv
import numpy
import tarfile
from importChemistryParameters import *



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


def getParam(paramName,file_conf):

    #leggo file conf
    # ##########################################
    #
    for line in open(file_conf,'r').readlines():
        thisline = line.split("=");

        if thisline[0]==paramName:
            return thisline[1].strip()


def importChemistry(lastGen,filesPath,ConfFilePath,outputFilePath):

    id_last_reaction=getLastReactionID(lastGen,filesPath)


    file_conf=ConfFilePath+'acsm2s.conf'
    lunghezza_max_fd=getParam('lunghezza_max_fd',file_conf)
    V1=float(getParam('volume',file_conf))
    InitialNumberOfSpecies=pow(2,int(lunghezza_max_fd)+1)-2

    overallConcentration=getParam('overallConcentration',file_conf)
    initialConc=(float(overallConcentration)/float(InitialNumberOfSpecies))*float(V1)

    #nome file spesis ultima generazione e ultima reazione
    file_species="%s%s%s%s%s%s" % (filesPath,'species_',lastGen,'_1_',id_last_reaction,'.csv')
    file_reactions="%s%s%s%s%s%s" % (filesPath,'reactions_',lastGen,'_1_',id_last_reaction,'.csv')
    file_catalysis="%s%s%s%s%s%s" % (filesPath,'catalysis_',lastGen,'_1_',id_last_reaction,'.csv')
    file_conf="%s%s" %(filesPath,'acsm2s.conf')
    file_influx="%s%s" %(filesPath,'_acsinflux.conf')

    file2write=outputFilePath+'ExportFile.dizzy'

    with open(file2write,'w') as file2write:

        volume= "V1 = "+ str(float(V1))+ ";\n"
        print volume
        file2write.write(volume)

        comment= "\n"+ "\n"+'//SPECIES'+ "\n"
        file2write.write(comment)

        print comment

        #costruisco dizionario con id specie e nome specie (e dizionario con nome specie e indicazione se bufferizzata) e scrivo elenco specie con concentrazioni su file
        speciesIDs={}
        speciesNames={}
        writtenComparments={} #siccome non posso specificare compartimento più di una volta ne tengo traccia

        for line in open(file_species,'r').readlines():
            thisline = line.split("\t");

            speciesID=thisline[0]
            speciesName=thisline[1]
            lockedConcentration=thisline[14]
            Kdecomp=str(float(thisline[5])/V1)


#            #se la concentrazione è bloccata la specie deve comparire nelle reazioni con il dollaro davanti
#            if (int(lockedConcentration) == 1):
#                speciesIDs[speciesID]='$'+speciesName
#            else:
#                speciesIDs[speciesID]=speciesName

            speciesIDs[speciesID]=speciesName
            speciesNames[speciesName]=lockedConcentration.strip()

            print 'speciesNames[speciesName] ',speciesName


            #se specie non era presente all'inzio metto concetrazione a 0
            if int(speciesID)<int(InitialNumberOfSpecies):
                initialConcentration=thisline[1]+'='+str(initialConc)+';'
            else:
                initialConcentration=thisline[1]+'='+str(0)+';'



            if thisline[1] not in writtenComparments:
                print initialConcentration
                file2write.write(initialConcentration)
                file2write.write('\n')

                speciesCompartment=thisline[1]+" @ V1;"+"\n"
                print speciesCompartment
                writtenComparments[thisline[1]]='V1'
                file2write.write(speciesCompartment)



        #costruisco dizionario con id reaction e specie catalyst(già tradotta dall'id) e rate
        catalysisIDs={}
        for line in open(file_catalysis,'r').readlines():
            thisline = line.split("\t");
            catalysisIDs[thisline[2]]= speciesIDs[thisline[1]]
            Kcleav=str(float(thisline[5])/V1)
            Kcond=str(float(thisline[4])/V1)
            Kcomp=str(float(thisline[6].strip())/V1)


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
                ##print 'substrate',substrate
                #se concentrazione è bloccata aggiungo dollaro
                ##print 'speciesNames[substrate]',speciesNames[substrate]
                if speciesNames[substrate]=='1':
                    ##print "==1"
                    substrate='$'+substrate

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
                    speciesCompartment=product+" @ V1;"+"\n"
                    file2write.write(speciesCompartment)
                    print speciesCompartment
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
                    speciesCompartment=thisline[1]+" @ V1;"+"\n"
                    file2write.write(speciesCompartment)
                    print speciesCompartment


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
