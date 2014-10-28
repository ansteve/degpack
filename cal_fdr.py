def cal_FDR(arr) :
        plen = len(arr)
        p_arr = np.copy(arr)
        q_arr = [1]*plen
        sort_ind = np.argsort(p_arr)
        p_arr.sort()
 
        for i in range(plen-1, -1, -1) :
                if (i==plen-1) :
                        q_arr[i]=p_arr[i]
                else :
                        rank = i+1
                        q_arr[i]=min(p_arr[i]*(float(plen)/float(rank)), q_arr[i+1])
        q = [0]*plen
        for i in range(0, plen) :
                q[sort_ind[i]] = q_arr[i]
