from collections import defaultdict
import numpy as np

def evalSpice(filename):
        # Initialize the dictionary to store components
        circuit_dict = defaultdict(list)

        # Initialize the nodes dictionary with GND
        nodes = {'GND': 0}
        in_circ=False
        try:
            open(filename,'r')
        except:
            raise FileNotFoundError('Please give the name of a valid SPICE file as input')

        with open(filename, 'r') as file:
            for line in file:
                words = line.strip().split()
                if not in_circ: #to dodge garbage before circuit
                    if ".circuit" in words:
                        in_circ=True
                    continue
                
                if ".end" in words: #to dodge garbage after the circuit
                    break
                        
                if len(words)<4: #if each line does not contain a minimum of 4 parameters
                    continue
                component_type = words[0][0] #extracting the first letter to check type of component
                #so my dictionary will look like:
                #{'voltage_sources':[{'name':v1,'node1:..},{'name:..}],'current_sources':[{'name':..},{}],..}
                if component_type == 'V':
                    for c in range(3,len(words)):
                        try:
                            float(words[c]) 
                            voltage=words[c]
                            break
                            
                            
                        except ValueError:
                            b1=0
                            
                    name, node1, node2 = words[0],words[1],words[2]
                    v1=voltage.split()
                    voltage=v1[len(v1)-1]
                    circuit_dict['voltage_sources'].append({'name': name, 'node1': node1, 'node2': node2, 'voltage': float(voltage)})
                elif component_type == 'I': 
                    for c in range(3,len(words)):
                        try:
                            float(words[c])
                            current=words[c]
                            break
                        except ValueError:
                            b1=0
                    name,node1,node2 = words[0],words[1],words[2]
                    circuit_dict['current_sources'].append({'name': name, 'node1': node1, 'node2': node2, 'current': float(current)})
                elif component_type == 'R':
                    name, node1, node2, resistance = words[0],words[1],words[2],words[3]
                    circuit_dict['resistances'].append({'name': name, 'node1': node1, 'node2': node2, 'resistance': float(resistance)})
                else:
                    raise ValueError('Only V, I, R elements are permitted')

                # Update the nodes dictionary
                nodes[node1] = None
                nodes[node2] = None
        if(not in_circ):
            raise ValueError('Malformed circuit file')

        # Assign consecutive numbers to nodes
        i = 1
        nodes['GND']=0
        for node in nodes.keys():
            if node != 'GND':
                nodes[node] = i
                i += 1
        #number of voltage sources        
        num_vs=len(circuit_dict['voltage_sources']) 
        #number of nodes
        num_nodes = len(nodes)
        num_is=len(circuit_dict['current_sources'])
        num_rs=len(circuit_dict['resistances'])
        if (num_is==0 and num_rs==0):
            raise ValueError('Circuit error: no solution')
        if(num_vs==0 and num_rs==0):
            raise ValueError('Circuit error: no solution')

        #creating the numpy array which will solve our equations, we need 2 matrices
        #the matrix dimentions will be num_nodes-1+num of voltage sources (beacuse we don't need ground node)
        dim=num_nodes+num_vs-1
        mat1=np.zeros((dim,1))
        mat2=np.zeros((dim,dim))
        #assigning values to matrix1
        #storing info about current flowing thru voltage sources in info_list
        info_list=[]
        for h in range(num_vs):
            info_list.append(0)
        #updating for finding values of current in branches with voltage sources
        for col in range(num_vs):
            nd1=nodes[circuit_dict['voltage_sources'][col]['node1']]
            nd2=nodes[circuit_dict['voltage_sources'][col]['node2']]
            if(nd1==0 or nd2==0):
                if(nd1==0):
                    mat2[nd2-1][col+num_nodes-1]-=1
                elif(nd2==0):
                    mat2[nd1-1][col+num_nodes-1]+=1
            else:
                mat2[nd1-1][col+num_nodes-1]+=1
                mat2[nd2-1][col+num_nodes-1]-=1
            info_list[col]=circuit_dict['voltage_sources'][col]['name']
           
        #adding voltage sources to matrix1 and handling the equation:difference in node potentials is equal to value of voltage connected across them    
        ctr = 0
        for j in range(num_nodes - 1, dim):
            mat1[j][0] = circuit_dict['voltage_sources'][ctr]['voltage']
            nd1=nodes[circuit_dict['voltage_sources'][ctr]['node1']]
            nd2=nodes[circuit_dict['voltage_sources'][ctr]['node2']]
            if(nd1==0 or nd2==0):
                if(nd1==0):
                    mat2[j][nd2-1]=-1
                    
                else:
                    mat2[j][nd1-1]=1
            else:
                mat2[j][nd1-1]=1
                mat2[j][nd2-1]=-1
            ctr=ctr+1
        #updating constant current sources in matrix1    
        for s1 in range(len(circuit_dict['current_sources'])):
            n1=nodes[circuit_dict['current_sources'][s1]['node1']]
            n2=nodes[circuit_dict['current_sources'][s1]['node2']]
            I_val=circuit_dict['current_sources'][s1]['current']
            if(n1==0 or n2==0):
                if(n1==0):
                    mat1[n2-1][0]+=I_val
                else:
                    mat1[n1-1][0]+=I_val  
            else:
                mat1[n1-1][0]+=I_val
                mat1[n2-1][0]-=I_val
        #populating the matrix with appropriate resistance values
        num_res=len(circuit_dict['resistances']) 
        for k in range(num_res):
            nod1=nodes[circuit_dict['resistances'][k]['node1']]
            nod2=nodes[circuit_dict['resistances'][k]['node2']]
            r_val=circuit_dict['resistances'][k]['resistance']
            if(nod1==0 or nod2==0):
                if(nod1==0):
                    mat2[nod2-1][nod2-1]+=(1/r_val)
                else:
                    mat2[nod1-1][nod1-1]+=(1/r_val)
            else:
                mat2[nod1-1][nod1-1]+=(1/r_val)
                mat2[nod2-1][nod2-1]+=(1/r_val)
                mat2[nod1-1][nod2-1]-=(1/r_val)
                mat2[nod2-1][nod1-1]-=(1/r_val)
        x=np.linalg.solve(mat2,mat1)
        dict1={} 
        #creating the dictionaries to be returned by appropriately de-indexing
        for w in nodes.keys():
            if w=='GND':
                dict1['GND']=0.0
            else:
                dict1[w]=x[nodes[w]-1][0]
        dict2={}
        cn=0
        for j in info_list:
            dict2[j]=x[cn+num_nodes-1][0]
            cn=cn+1
        
        
        
        return dict1,dict2


  
