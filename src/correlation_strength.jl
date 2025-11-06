
module correlation_strength

using JuMP, DataFrames, LinearAlgebra, Leiden

export orcadot, NLPCorrStrengGenerating
function reduced(A,red)
    for i in 1:size(A,1)
        for j in 1:size(A,1)
            A[i,j]=round(A[i,j]*1000,digits=0) # make the numbers be round in the first digit
        end
    end
    #println(A)
    global res=1
    global ov_ccom=Leiden.leiden(A,resolution=res)
    global ov_cpart=ov_ccom[2] # the partition vector
    #run leiden until resolution gives no communities
    while length(ov_cpart)!=size(A,1)
        global ov_ccom=Leiden.leiden(A,resolution=res)
        global ov_cpart=ov_ccom[2]
        global res+=1 #when all the elements are divided individually, loop stops.
    end   # but at this time the value of res is actually one larger than the res that all elements are divided individually
    #run leiden until "red" number of objective communities is reached
    while ( length(ov_cpart)>red && res>0 ) #
        global ov_ccom=Leiden.leiden(A,resolution=res)
        global ov_cpart=ov_ccom[2]
        global res-=.1
    end

    return(ov_cpart) # just give out a value that has less dimension than red, but I think that just should be the integer lower than red, can't figure out the meaning of it...
end

function NLPCorrStrengGenerating(IneqMatrix,ObjMatrix,EqMatrix,NoParitition,NoOfSP=1) # version 8/2/2024 consider the condition that EqMatrix is all zero
    # IneqMatrix, Jacobian matrix of inequal constaints; ObjMatrix, Jacobian matrix of objectives; EqMatrix, Jacobian matrix of equal constaints; NoParitition and NoOfSP, No. of paritions and selected points
    NoOfObj = Int(size(ObjMatrix,1)/NoOfSP) # number of selected points times number of obectives equals the numbers objectives in the ObjMatrix
    # println("NoOfObj = ", NoOfObj)
    normc=zeros(size(ObjMatrix)) # normalized Cmat matrix
    for i in 1:size(ObjMatrix,1)
      normc[i,:]=-normalize(ObjMatrix[i,:])#normalize the value of the matrix row by row
    end # thus normc is the normalized ce
    posmin=0.9 # hyperparameters in weight function usually 0.9
    beta = 100 # hyperparameters in weight function usually 100
    IneqSecstrength=zeros(NoOfObj,NoOfObj,size(IneqMatrix,1),NoOfSP) #Sijkn here for inequal constraints
    totIneqWt=zeros(NoOfObj,NoOfObj,size(IneqMatrix,1),NoOfSP) #Wijkn here for inequal constraints
    #Inequality constraints
    for i in 1:NoOfObj
      for j in i+1:NoOfObj
        for k in 1:size(IneqMatrix,1)
          for n in 1:NoOfSP
            # println("principle one: ", dot(normc[i,:],Ae[k,:]))
            # println("principle two: ", dot(normc[j,:],Ae[k,:]))
            cp1=normc[(i-1)*NoOfSP+n,:]-dot(normc[(i-1)*NoOfSP+n,:],IneqMatrix[k,:])/(sqrt(sum(IneqMatrix[k,p]^2 for p in 1:size(IneqMatrix,2)))^2)*IneqMatrix[k,:] # to project normc[i] on Ae[K]
            cn1=normc[(i-1)*NoOfSP+n,:]-cp1 # get the left part that is normal to Ae[k]
            cp2=normc[(j-1)*NoOfSP+n,:]-dot(normc[(j-1)*NoOfSP+n,:],IneqMatrix[k,:])/(sqrt(sum(IneqMatrix[k,p]^2 for p in 1:size(IneqMatrix,2)))^2)*IneqMatrix[k,:]
            cn2=normc[(j-1)*NoOfSP+n,:]-cp2
            if cp1!=zeros(length(cp1)) 
                cp1norm=normalize(cp1) # normalize the proejcted vectors
            else
                cp1norm=cp1 # function narmalize can be applied to all zeros!!!
            end
            if cp2!=zeros(length(cp2))
                cp2norm=normalize(cp2)
            else
                cp2norm=cp2
            end
            if dot(normc[(i-1)*NoOfSP+n,:],IneqMatrix[k,:])!=0 && dot(normc[(j-1)*NoOfSP+n,:],IneqMatrix[k,:])!=0 # if this two objective vectors are not normal to a row in the constrain matrix
              #Previously we don't count the objecives that are pointing into the contraints because we can only have optimal solutions in the boundary of the feasible region in linear problem.
              #But now we are considering nonlinear problem, we need to consider all the objectives. Optimal points can be in the middle of the feasible region when objectives are nonliear.
              if dot(IneqMatrix[k,:],cn1)>0 && dot(IneqMatrix[k,:],cn2)>0 # normally the dot product results should be zero because the cni is normal to Ae[k]
                IneqSecstrength[i,j,k,n] = dot(cp1norm,cp2norm)
                totIneqWt[i,j,k,n] = 1-posmin*(1/(1+exp(-beta*IneqSecstrength[i,j,k,n]))) # calculat the weight Sij*Wijk
                # if secstrength==0
                #   totIneqWt[i,j,k]=0 # why
                # else
                #   totIneqWt[i,j,k]+=(1-posmin*(1/(1+exp(-beta*secstrength)))) # store the Wijk
                # end
              else
                # println("one here")siz
                IneqSecstrength[i,j,k,n]=0
                totIneqWt[i,j,k,n] = 0
              end
            else#if one objective is perellel to the constraint, I think we should still calculate it becuase its also in the boundary and there should be tradeoffs.
              IneqSecstrength[i,j,k,n] = dot(cp1norm,cp2norm)
              totIneqWt[i,j,k,n] = 1-posmin*(1/(1+exp(-beta*IneqSecstrength[i,j,k,n])))
            end
            #println(cp1,cp2)
            #println(cp1norm, cp2norm)
            #println(secstrength, " secondary strength")
            #dsm[i,j]+=(1-(secstrength+1)*0.5)*secstrength#secstrength
          end
        end
      end
    end
    ##############################################################################
    #Equality constraints
    EqSecstrength=zeros(NoOfObj,NoOfObj,size(EqMatrix,1),NoOfSP) #Sijkn here for equal constraints
    totEqWt=zeros(NoOfObj,NoOfObj,size(EqMatrix,1),NoOfSP) #Wijkn here for equal constraints
    for i in 1:NoOfObj
      for j in i+1:NoOfObj
        for k in 1:size(EqMatrix,1)
          for n in 1:NoOfSP
            if sum(EqMatrix) !=0 # if there is no equal constraint in the problem, we just make the weights are zero. 
              #Also there will be NaN in the cps if this if statement is not added..
              cp1=normc[(i-1)*NoOfSP+n,:]-dot(normc[(i-1)*NoOfSP+n,:],EqMatrix[k,:])/(sqrt(sum(EqMatrix[k,p]^2 for p in 1:size(EqMatrix,2)))^2)*EqMatrix[k,:]
              cp2=normc[(j-1)*NoOfSP+n,:]-dot(normc[(j-1)*NoOfSP+n,:],EqMatrix[k,:])/(sqrt(sum(EqMatrix[k,p]^2 for p in 1:size(EqMatrix,2)))^2)*EqMatrix[k,:]
      
              if cp1!=zeros(length(cp1))
                  cp1norm=normalize(cp1)
              else
                  cp1norm=cp1
              end
              if cp2!=zeros(length(cp2))
                  cp2norm=normalize(cp2)
              else
                  cp2norm=cp2
              end
              if dot(normc[(i-1)*NoOfSP+n,:],EqMatrix[k,:])!=0 && dot(normc[(j-1)*NoOfSP+n,:],EqMatrix[k,:])!=0
                EqSecstrength[i,j,k,n] = dot(cp1norm,cp2norm)
                totEqWt[i,j,k,n] = 1
                #   if seceqstren==0
                #       toteqwt[i,j,k]=0
                #   else
                #       toteqwt[i,j,k]+=(1-posmin*(1/(1+exp(-beta*seceqstren)))) # might be a bug
                #   end
                # println("Wijk: ", toteqwt[i,j,k]) # No abnorm
              else
                EqSecstrength[i,j,k,n] = dot(cp1norm,cp2norm)
                totEqWt[i,j,k,n] = 1
              end
            else
              EqSecstrength[i,j,k,n] = 0
              totEqWt[i,j,k,n] = 0
            end
          end
        end
      end
    end
    totwt= zeros(NoOfObj,NoOfObj)
    totSecWt=zeros(NoOfObj,NoOfObj)
    #deltaine=ones(size(ce,1),size(ce,1),size(Ae,1))
    #deltaeq=ones(size(ce,1),size(ce,1),size(de,1))
    for i in 1:NoOfObj
      for j in i+1:NoOfObj
        totwt[i,j]=totwt[j,i]= sum(totIneqWt[i,j,:,:])+sum(totEqWt[i,j,:,:])#sum of the total weight 
        totSecWt[i,j]=totSecWt[j,i]+=sum(totIneqWt[i,j,:,:].*IneqSecstrength[i,j,:,:])+sum(totEqWt[i,j,:,:].*EqSecstrength[i,j,:,:])
      end
    end
    adjMatrix=zeros(NoOfObj,NoOfObj)
    combmin=zeros(NoOfObj,NoOfObj)
    for i in 1:NoOfObj
      for j in i+1:NoOfObj
        adjMatrix[i,j]=adjMatrix[j,i]=(1/2)*(1 + totSecWt[i,j]/totwt[i,j])
        # combmin[i,j]=combmin[j,i]=minimum(vcat(dsmeq[i,j,:],dsm[i,j,:])) #debuging line.
      end
    end
    # for i in 1:NoOfObj
    #     for j in i+1:NoOfObj
    #         if isnan(comb[i,j])
    #             comb[i,j]=comb[j,i]=0
    #         end
    #     end
    # end
    # for i in 1:NoOfObj
    #     for j in i+1:NoOfObj
    #         comb[i,j]=comb[j,i]=(1/2)*(comb[i,j]+1)# result of Sij
    #         # combmin[i,j]=combmin[j,i]=(1/2)*(combmin[i,j]+1)
    #     end
    # end
    for i in 1:NoOfObj
      adjMatrix[i,i] = 1
    end
    println(adjMatrix)
    groups=reduced(copy(adjMatrix),NoParitition)
    
    #dpm, dsmfin, dsmeqfin, groups,
    return (adjMatrix = adjMatrix, totwt = totwt, groups = groups, totIneqWt = totIneqWt)#, dsmfin, dsmeqfin, dpm
end
function orcadot(Ae,ce,de,comm)# Ae is the constrain matrix and ce is the objective matrix
  normc=zeros(size(ce)) # identical size of ce but fill with all 0 values.
  posmin=1
  for i in 1:size(ce,1)
    normc[i,:]=-normalize(ce[i,:])#normalize the value of the matrix row by row
  end # thus normc is the normalized ce
  #println(normc)
  dpm=zeros(size(ce,1),size(ce,1)) # might be different if ce is not a square. Thus, make a matrix of 0 values with same # of rows and columns of the # of ce's rows
  for i in 1:size(ce,1)
    for j in i+1:size(ce,1) # the # of columes at least should be larger than the # of rows.
      dpm[i,j]=dot(normc[i,:],normc[j,:]) # the upper triangle is filled with the dots products of rows in normc, of which the length of the vectors is 1. 
    end # all the dot product results of different rows in normc are stored in dmp.
  end
  dsm=zeros(size(ce,1),size(ce,1),size(Ae,1))
  totineqwt=zeros(size(ce,1),size(ce,1),size(Ae,1))
  deltaine=zeros(size(ce,1),size(ce,1),size(Ae,1))
  #Inequality constraints
  for i in 1:size(ce,1)
    for j in i+1:size(ce,1)
      for k in 1:size(Ae,1)
    if dot(normc[i,:],Ae[k,:])!=0 && dot(normc[j,:],Ae[k,:])!=0 # if this two objective vectors are not normal to a row in the constrain matrix
        cp1=normc[i,:]-dot(normc[i,:],Ae[k,:])/(sqrt(sum(Ae[k,p]^2 for p in 1:size(Ae,2)))^2)*Ae[k,:] # to project normc[i] on Ae[K]
        cn1=normc[i,:]-cp1 # get the left part that is normal to Ae[k]
        cp2=normc[j,:]-dot(normc[j,:],Ae[k,:])/(sqrt(sum(Ae[k,p]^2 for p in 1:size(Ae,2)))^2)*Ae[k,:]
        cn2=normc[j,:]-cp2
        if cp1!=zeros(length(cp1)) 
            cp1norm=normalize(cp1) # normalize the proejcted vectors
        else
            cp1norm=cp1 # mean if cp1 is not all zero?, which means that normc[i] is acutually normal to Ae[k]
        end
        if cp2!=zeros(length(cp2))
            cp2norm=normalize(cp2)
        else
            cp2norm=cp2
        end

        if dot(Ae[k,:],cn1)>0 && dot(Ae[k,:],cn2)>0 # normally the dot product results should be zero because the cni is normal to Ae[k]
            secstrength=dot(cp1norm,cp2norm)
            deltaine[i,j,k]=secstrength#abs(secstrength-dpm[i,j])  #stroage of Sij
            dsm[i,j,k]+=secstrength*(1-posmin*(1/(1+exp(-10*secstrength)))) # calculat the weight Sij*Wijk
            if secstrength==0
                totineqwt[i,j,k]=0 # why
            else
                totineqwt[i,j,k]+=(1-posmin*(1/(1+exp(-10*secstrength)))) # store the Wijk
            end
        else
            secstrength=0
            dsm[i,j,k]+= secstrength
            totineqwt[i,j,k]+=0
        end
    else
         secstrength=0
        dsm[i,j,k]+= secstrength
        totineqwt[i,j,k]+=0
    end
        #println(cp1,cp2)
        #println(cp1norm, cp2norm)

        #println(secstrength, " secondary strength")
        #dsm[i,j]+=(1-(secstrength+1)*0.5)*secstrength#secstrength

      end
    end
  end
  deltaine=deltaine/maximum(deltaine) # rescale the biggest as one..
  ##############################################################################
  #Equality constraints
  dsmeq=zeros(size(ce,1),size(ce,1),size(de,1))
  toteqwt=zeros(size(ce,1),size(ce,1),size(de,1))
  deltaeq=zeros(size(ce,1),size(ce,1),size(de,1))
  for i in 1:size(ce,1)
    for j in i+1:size(ce,1)
      for k in 1:size(de,1)
        cp1=normc[i,:]-dot(normc[i,:],de[k,:])/(sqrt(sum(de[k,p]^2 for p in 1:size(de,2)))^2)*de[k,:]
        cp2=normc[j,:]-dot(normc[j,:],de[k,:])/(sqrt(sum(de[k,p]^2 for p in 1:size(de,2)))^2)*de[k,:]

        if cp1!=zeros(length(cp1))
            cp1norm=normalize(cp1)
        else
            cp1norm=cp1
        end
        if cp2!=zeros(length(cp2))
            cp2norm=normalize(cp2)
        else
            cp2norm=cp2
        end
        if dot(normc[i,:],de[k,:])!=0 && dot(normc[j,:],de[k,:])!=0

            seceqstren=dot(cp1norm,cp2norm)
            deltaeq[k]=abs(seceqstren-dpm[i,j])
            dsmeq[i,j,k]+=seceqstren*1#(1-posmin*(1/(1+exp(-10*seceqstren))))
            if seceqstren==0
                toteqwt[i,j,k]=0
            else
                toteqwt[i,j,k]+=(1-posmin*(1/(1+exp(-10*seceqstren)))) # might be a bug
            end
        else
            seceqstren=0
            deltaeq[i,j,k]=0
            dsmeq[i,j,k]+=0
            toteqwt[i,j,k]+=0
        end
      end
    end
  end
  deltaeq=deltaeq/maximum(deltaeq)


  totwt= zeros(size(ce,1),size(ce,1))
  totdsm=zeros(size(ce,1),size(ce,1))
  #deltaine=ones(size(ce,1),size(ce,1),size(Ae,1))
  #deltaeq=ones(size(ce,1),size(ce,1),size(de,1))
  for i in 1:size(ce,1)
    for j in i+1:size(ce,1)
      totwt[i,j]=totwt[j,i]+=sum(totineqwt[i,j,k]*deltaine[i,j,k] for k in 1:size(Ae,1))+sum(toteqwt[i,j,k]*deltaeq[i,j,k] for k in 1:size(de,1))# upper part of Sij
      totdsm[i,j]=totdsm[j,i]+=sum(dsmeq[i,j,k]*deltaeq[i,j,k] for k in 1:size(de,1))+sum(dsm[i,j,k]*deltaine[i,j,k] for k in 1:size(Ae,1))
      #if dsm[i,j] ==0 || totineqwt[i,j]==0
    #    totineqwt[i,j]=1
     # end
      #if dsmeq[i,j]==0 || toteqwtz[i,j]==0
    #    toteqwt[i,j]=1
     # end
    end
  end
  #println(dsm, " ", totineqwt)
  #println(dsmeq," ", toteqwt)
  #dsmfin=zeros(size(ce,1),size(ce,1))
  #dsmeqfin=zeros(size(ce,1),size(ce,1))
  comb=zeros(size(ce,1),size(ce,1))
  combmin=zeros(size(ce,1),size(ce,1))
  for i in 1:size(ce,1)
    for j in i+1:size(ce,1)
      #dsmfin[i,j]=dsm[i,j]/totineqwt[i,j]
      #dsmeqfin[i,j]=dsmeq[i,j]/toteqwt[i,j]
      comb[i,j]=comb[j,i]=totdsm[i,j]/totwt[i,j]# why is like this
      combmin[i,j]=combmin[j,i]=minimum(vcat(dsmeq[i,j,:],dsm[i,j,:])) #debuging line.
    end
  end
  for i in 1:size(ce,1)
      for j in i+1:size(ce,1)
          if isnan(comb[i,j])
              comb[i,j]=comb[j,i]=0
          end
      end
  end

  for i in 1:size(ce,1)
      for j in i+1:size(ce,1)
          comb[i,j]=comb[j,i]=(1/2)*(comb[i,j]+1)# result of Sij
          combmin[i,j]=combmin[j,i]=(1/2)*(combmin[i,j]+1)

      end
  end
  groups=reduced(copy(comb),comm)
  #dpm, dsmfin, dsmeqfin, groups,
  return comb,  combmin, groups, dsmeq, dsm, dpm, deltaeq, deltaine#, dsmfin, dsmeqfin, dpm
end

end # module correlation_strength