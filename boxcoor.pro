function boxcoor,centerX,centerY,sizeX=sizeX,sizeY=sizeY,fordraw=fordraw,wavelength=wavelength
;; calculates the corners of the box given the center points
;; returns the corners (x1,x2,y1,y2)
;; fordraw returns x,y points in order to draw boxes

x1 = double(centerX) - double(sizeX)/2d
x2 = double(centerX) + double(sizeX)/2d
y1 = double(centerY) - double(sizeY)/2d
y2 = double(centerY) + double(sizeY)/2d

if keyword_set(fordraw) then begin
   ;; half all points for the half-size image
   x1 /= 2d
   x2 /= 2d
   y1 /= 2d
   y2 /= 2d
   return,[[x1,y1],[x2,y1],[x2,y2],[x1,y2],[x1,y1]]
endif else begin
   if keyword_set(wavelength) then begin
      assert,centerY,'>=',375,"Extraction box is below the acceptable range!"
      assert,centerY,'<=',1015,"Extraction box is above the acceptable range!"

      if centerY GT 880 then begin ;; order 3
         lambdaL = 2.46D ;; in microns
         lambdaR = 1.88D ;; in microns
      endif else begin
         if centerY GT 725 and centerY LE 880 then begin ;; order 4
            lambdaL = 1.85D
            lambdaR = 1.42D ;; in microns
         endif else begin
            if centerY GT 590 and centerY LE 725 then begin ;; order 5
               lambdaL = 1.48D
               lambdaR = 1.13D
            endif else begin
               if centerY GT 375 and centerY LE 590 then begin ;; order 6
                  lambdaL = 1.23D
                  lambdaR = 0.95D
               endif
            endelse
         endelse
      endelse
            
      ;;given the x coordinate in pixels, find the wavelength in
      ;;microns, set by the points (2048,lambdaR) and (0,lambdaL)
      lam = lambdaL + (centerX - 2048d) * (lambdaR - lambdaL)/(0d - 2048d)
      return,lam
   endif else return, [[x1,y1],[x2,y2]]
endelse

end
