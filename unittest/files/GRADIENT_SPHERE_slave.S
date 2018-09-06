        .set noat
        .set noreorder
        .arch sw3
        .section .rodata, "a", "progbits"
        .align 4
        .section .text1, "ax", "progbits"
        .align 4
        .align 4
        .globl GRADIENT_SPHERE 
        .ent GRADIENT_SPHERE
GRADIENT_SPHERE:	     
             vshff $9,$10,0x00,$9          
		     vshff $11,$12,0x00,$11        
             vshff $13,$14,0x00,$13        
		     vshff $15,$26,0x00,$26     
             vshff $9,$11,0x88,$28         
	         vshff $13,$26,0x88,$26        		   
		     #LD IADDR1 ARRAY           
             vldd $0,0($17)
			 vldd $1,32($17)                 
             vldd $2,64($17)                 
             vldd $3,96($17) 	          
             #LD IADDR2 ARRAY
			 vldd $4,0($18)               	
             vldd $5,32($18) 	          
	 	     vldd $6,64($18)               
		     vldd $7,96($18)
 		
             #vshff IADDR1 ARRAY
             #vshff $1,$0,0x44,$8
			 #vshff $1,$0,0xee,$9
			 #vshff $3,$2,0x44,$10
			 #vshff $3,$2,0xee,$11
			 #vshff $10,$8,0x88,$12
			 #vshff $10,$8,0xdd,$13
			 #vshff $11,$9,0x88,$14
			 #vshff $11,$9,0xdd,$15
			 vshff $4,$4,0x0,$17
			 vshff $4,$4,0x55,$18
			 vshff $4,$4,0xaa,$21
			 vshff $4,$4,0xff,$22
             vshff $5,$5,0x0,$23
			 vshff $5,$5,0x55,$24
			 vshff $5,$5,0xaa,$25
		     vshff $5,$5,0xff,$27

			 vmuld $17,$0,$8
			 vmuld $23,$0,$9
			 vmad  $18,$1,$8,$8
			 vmad  $24,$1,$9,$9
			 vmad  $21,$2,$8,$8
			 vmad  $25,$2,$9,$9
			 vmad  $22,$3,$8,$8
			 vmad  $27,$3,$9,$9
			 vshff $6,$6,0x0,$17
			 vshff $6,$6,0x55,$18
		       vshff $6,$6,0xaa,$21
		       vshff $6,$6,0xff,$22
		       vshff $7,$7,0x0,$23
		       vshff $7,$7,0x55,$24
		       vshff $7,$7,0xaa,$25
		       vshff $7,$7,0xff,$27
			 vmuld $17,$0,$10
			 vmuld $23,$0,$11
			 vmad  $18,$1,$10,$10
			 vmad  $24,$1,$11,$11
			 vmad  $21,$2,$10,$10
			 vmad  $25,$2,$11,$11
			 vmad  $22,$3,$10,$10
			 vmad  $27,$3,$11,$11
			 
			 vshff $20,$20,0x0,$20
  
			 #vshff IADDR2 ARRAY
			 vshff $5,$4,0x44,$17
			 vshff $5,$4,0xee,$18
			 vshff $7,$6,0x44,$21
			 vshff $7,$6,0xee,$22
			 vshff $21,$17,0x88,$12
			 vshff $21,$17,0xdd,$13
			 vshff $22,$18,0x88,$14
			 vshff $22,$18,0xdd,$15
             
             vshff $0,$0,0x0,$23
             vshff $0,$0,0x55,$24
             vshff $0,$0,0xaa,$25
             vshff $0,$0,0xff,$27
             vshff $1,$1,0x0,$4
             vshff $1,$1,0x55,$5
             vshff $1,$1,0xaa,$6
             vshff $1,$1,0xff,$7
             vmuld $23,$12,$21
             vmuld $4,$12,$22
             vmad  $24,$13,$21,$21
             vmad  $5,$13,$22,$22
             vmad  $25,$14,$21,$21
             vmad  $6,$14,$22,$22
             vmad  $27,$15,$21,$21
             vmad  $7,$15,$22,$22

            
             vshff $2,$2,0x0,$0
             vshff $2,$2,0x55,$1
             vshff $2,$2,0xaa,$17
             vshff $2,$2,0xff,$18
             vshff $3,$3,0x0,$4
             vshff $3,$3,0x55,$5
             vshff $3,$3,0xaa,$6
             vshff $3,$3,0xff,$7
             vmuld $0,$12,$23
             vmuld $4,$12,$24
             vmad  $1,$13,$23,$23
             vmad  $5,$13,$24,$24
             vmad  $17,$14,$23,$23
             vmad  $6,$14,$24,$24
             vmad  $18,$15,$23,$23
             vmad  $7,$15,$24,$24
		 
		 vldd $0,0($19)
             vldd $1,32($19)
             vldd $2,64($19)
             vldd $3,96($19)
             vldd $12,128($19)
             vldd $13,160($19)
             vldd $14,192($19)
             vldd $15,224($19)
            
             vldd $4,256($19)
             vldd $5,288($19)
             vldd $6,320($19)
             vldd $7,352($19)
             vldd $17,384($19)
             vldd $18,416($19)
             vldd $25,448($19)
             vldd $27,480($19)
             
             vmuld $0,$21,$0
             vmuld $1,$22,$1
             vmuld $2,$23,$2
		 vmuld $3,$24,$3
             vmuld $4,$21,$4
             vmuld $5,$22,$5
             vmuld $6,$23,$6
             vmuld $7,$24,$7 
           

           vmad  $12,$8,$0,$0
           vmad  $13,$9,$1,$1
           vmad  $14,$10,$2,$2
           vmad  $15,$11,$3,$3
		vmad  $17,$8,$4,$4
		vmad  $18,$9,$5,$5
		vmad  $25,$10,$6,$6
		vmad  $27,$11,$7,$7
			
		 vmuld $0,$20,$0
             vmuld $1,$20,$1
             vmuld $2,$20,$2
             vmuld $3,$20,$3
             vmuld $4,$20,$4
             vmuld $5,$20,$5
             vmuld $6,$20,$6
		 vmuld $7,$20,$7
		
               
		 vstd  $0,0($16)              
		 vstd  $1,32($16)             
		 vstd  $2,64($16)             
		 vstd  $3,96($16)
             vstd  $4,128($16)
             vstd  $5,160($16)
             vstd  $6,192($16)
             vstd  $7,224($16)             
             vextf $27,0x0,$9               
             vextf $27,0x1,$10              
		 vextf $27,0x2,$11              
             vextf $27,0x3,$12              
             vextf $26,0x1,$13              
             vextf $26,0x2,$14              
		 vextf $26,0x3,$15    
             ret  
                     .end GRADIENT_SPHERE
#ifdef __slave_alias__
        .weak slave_GRADIENT_SPHERE;slave_GRADIENT_SPHERE=GRADIENT_SPHERE
        .weak slave_gradient_sphere;slave_gradient_sphere=GRADIENT_SPHERE
		.weak slave_gradient_sphere_;slave_gradient_sphere_=GRADIENT_SPHERE
#endif
         
