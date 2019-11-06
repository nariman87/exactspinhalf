function Y = CycleBits(state, N)
 
 State = de2bi(state);
 State = fliplr(State);
 
 State_N = zeros(1, N);
 State_rest = zeros( 1 , max(0, size(State,2) - N) );
 
 for i = 1 : size(State,2)
    if ( i > size(State,2) - N )
       State_N( i - (size(State,2) - N) ) = State(i);
    else
       State_rest(i) = State(i);      
    end
 end
 
 State_N = ( circshift(State_N',-1) )';
 State = [ State_rest State_N ];
  
 State = fliplr(State);
 Y = bi2de(State);
 
end
