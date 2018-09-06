        .set noat
        .set noreorder
        .arch sw3
        .section .rodata, "a", "progbits"
        .align 4
        .section .text1, "ax", "progbits"
        .align 4
        .align 4
        .globl VORTICITY_SPHERE 
        .ent VORTICITY_SPHERE
VORTICITY_SPHERE:	     
                     vshff $9,$10,0x00,$9          
		     vshff $11,$12,0x00,$11        
                     vshff $13,$14,0x00,$13        
		     vshff $15,$28,0x00,$28        
                     vshff $9,$11,0x88,$27         
	             vshff $13,$28,0x88,$28        		   
		     #LD IADDR2 ARRAY           
                     vldd $0,0($18)                  
                     vldd $1,32($18)                 
                     vldd $2,64($18)                 
                     vldd $3,96($18) 	          
                     vldd $4,128($18)               	
                     vldd $5,160($18) 	          
	 	     vldd $6,192($18)               
		     vldd $7,224($18)

		     #LD IADDR1(1) ARRAY          
		     vldd $8,0($17)                 
		     vldd $9,32($17)                
		     vldd $10,64($17)                
                     vldd $11,96($17)    		  
                     vldd $12,128($17)  		  
                     vldd $13,160($17)               
                     vldd $14,192($17)               
                     vldd $15,224($17)  
                     
                     vmuld $0,$8,$8
                     vmuld $1,$9,$9
                     vmuld $2,$10,$10
                     vmuld $3,$11,$11
		     vmad $4,$12,$8,$8
                     vmad $5,$13,$9,$9
                     vmad $6,$14,$10,$10
                     vmad $7,$15,$11,$11
 


                     vldd $22,256($17) 
                     vldd $23,288($17) 
                     vldd $24,320($17)     
                     vldd $25,352($17)
                     vldd $12,384($17)
                     vldd $13,416($17)
                     vldd $14,448($17)
                     vldd $15,480($17)
	          
                     #VCO(2)=IADDR1(1)*IADDR2(1)+IADDR1(2)*IADDR2(2)    		
                     vmuld $0,$22,$0                
                     vmuld $1,$23,$1   	          
                     vmuld $2,$24,$2   		  
                     vmuld $3,$25,$3                
		     vmad $4,$12,$0,$4             
		     vmad $5,$13,$1,$5             
		     vmad $6,$14,$2,$6             
		     vmad $7,$15,$3,$7 
		      
		     #LD IADDR3 ARRAY               
		     vldd $22,0($19)                 
		     vldd $23,32($19)                
		     vldd $24,64($19)                
		     vldd $25,96($19)                
		     #vshff IADDR3                  
		     vshff $22,$22,0x0,$0          
		     vshff $22,$22,0x55,$1          
		     vshff $22,$22,0xaa,$2          
		     vshff $22,$22,0xff,$3          
		     vshff $23,$23,0x0,$12          
		     vshff $23,$23,0x55,$13          
		     vshff $23,$23,0xaa,$14         
		     vshff $23,$23,0xff,$15          
		     #VTEMP=VCO*IADDR3              
		     vmuld $0,$8,$0               
		     vmuld $12,$8,$12                   
		     vmad  $1,$9,$0,$0               
		     vmad  $13,$9,$12,$12               
 		     vmad  $2,$10,$0,$0                 
		     vmad  $14,$10,$12,$12                
		     vmad  $3,$11,$0,$0                  
		     vmad  $15,$11,$12,$12
		     vshff $24,$24,0x0,$17         
		     vshff $24,$24,0x55,$18         
		     vshff $24,$24,0xaa,$19        
		     vshff $24,$24,0xff,$3        
		     vshff $25,$25,0x0,$2        
		     vshff $25,$25,0x55,$13        
		     vshff $25,$25,0xaa,$14        
		     vshff $25,$25,0xff,$15       
		     vmuld $17,$8,$1             
		     vmuld $2,$8,$2             
		     vmad  $18,$9,$1,$1           
		     vmad  $13,$9,$2,$2          
 		     vmad  $19,$10,$1,$1        
		     vmad  $14,$10,$2,$2          
		     vmad  $3,$11,$1,$1         
		     vmad  $15,$11,$2,$2 
                     #VSHFF IADDR3                  
		     vshff $23,$22,0x44,$17       
		     vshff $23,$22,0xee,$18       
		     vshff $25,$24,0x44,$13       
                     vshff $25,$24,0xee,$19       
                     vshff $13,$17,0x88,$8         
                     vshff $13,$17,0xdd,$9         
                     vshff $19,$18,0x88,$10         
                     vshff $19,$18,0xdd,$11        
                     #vshff VTEMP               
                     vshff $4,$4,0x0,$22  	  
                     vshff $4,$4,0x55,$23        
                     vshff $4,$4,0xaa,$24     	  
                     vshff $4,$4,0xff,$25        
                     vshff $5,$5,0x0,$13   	  
                     vshff $5,$5,0x55,$14        
                     vshff $5,$5,0xaa,$15     	  
                     vshff $5,$5,0xff,$17        
                     vmuld $22,$8,$22           
                     vmuld $13,$8,$13            
                     vmad  $23,$9,$22,$22        
                     vmad  $14,$9,$13,$13         
                     vmad  $24,$10,$22,$22          
                     vmad  $15,$10,$13,$13          
                     vmad  $25,$11,$22,$22         
                     vmad  $17,$11,$13,$13          
                     
		     vsubd $22,$0,$4               
		     vsubd $13,$12,$5               
		     vshff $6,$6,0x0,$25   	  
                     vshff $6,$6,0x55,$24        
                     vshff $6,$6,0xaa,$23     	  
                     vshff $6,$6,0xff,$22        
                     vshff $7,$7,0x0,$15   	  
                     vshff $7,$7,0x55,$14        
                     vshff $7,$7,0xaa,$13     	  
                     vshff $7,$7,0xff,$12        
                     vmuld $25,$8,$25          
                     vmuld $15,$8,$15            
                     vmad  $24,$9,$25,$25         
                     vmad  $14,$9,$15,$15         
                     vmad  $23,$10,$25,$25          
                     vmad  $13,$10,$15,$15          
                     vmad  $22,$11,$25,$25         
                     vmad  $12,$11,$15,$15

                     vsubd $25,$1,$6               
		     vsubd $15,$2,$7
                             
		     vshff $21,$21,0x0,$21              
		     vldd $8,0($20)                 
		     vldd $9,32($20)          
		     vldd $10,64($20)          
		     vldd $11,96($20)          
		     vmuld $21,$4,$4          
		     vmuld $21,$5,$5          
		     vmuld $21,$6,$6        
		     vmuld $21,$7,$7        
		     vmuld $8,$4,$4        
		     vmuld $9,$5,$5        
		     vmuld $10,$6,$6        
		     vmuld $11,$7,$7        
		     vstd  $4,0($16)              
		     vstd  $5,32($16)             
		     vstd  $6,64($16)             
		     vstd  $7,96($16)             
                     vextf $27,0x0,$9               
                     vextf $27,0x1,$10              
		     vextf $27,0x2,$11              
                     vextf $27,0x3,$12              
                     vextf $28,0x1,$13              
                     vextf $28,0x2,$14              
		     vextf $28,0x3,$15    
                     ret  
                     .end VORTICITY_SPHERE
#ifdef __slave_alias__
        .weak slave_VORTICITY_SPHERE;slave_VORTICITY_SPHERE=VORTICITY_SPHERE
        .weak slave_vorticity_sphere;slave_vorticity_sphere=VORTICITY_SPHERE
		.weak slave_vorticity_sphere_;slave_vorticity_sphere_=VORTICITY_SPHERE
		.weak slave_DIVERGENCE_SPHERE;slave_DIVERGENCE_SPHERE=VORTICITY_SPHERE
		.weak slave_divergence_sphere;slave_divergence_sphere=VORTICITY_SPHERE
		.weak slave_divergence_sphere_;slave_divergence_sphere_=VORTICITY_SPHERE
#endif
         
