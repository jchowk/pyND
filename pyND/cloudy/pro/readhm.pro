PRO readhm
  output_filename = 'hmr_uvb.fits'

  openr, unit, 'hmr_uvb_all.dat', /get_lun

  input_string = ''
  readf, unit, input_string
  zstring = str_sep(strtrim(strcompress(input_string), 2), ' ')
  z = float(zstring[2:n_elements(zstring)-1])

  num_redshifts = n_elements(z)

  hmspec = {z:0., wave:fltarr(167), energy:fltarr(167), $
            fnu:fltarr(167), comments:''}
  hmspec = replicate(hmspec, num_redshifts)
  
  hmspec.z = z
  hmspec.comments = 'Wavelength [Angstrom]; Energy [Rydbergs]; ' + $
     'F_nu [ergs/cm^2/s/Hz/sr]'

  i = 0
  readf, unit, input_string
  WHILE strmid(input_string, 0, 3) NE 'EOF' DO BEGIN
    input_string_array = $
       str_sep(strtrim(strcompress(input_string), 2), ' ')
    num_entries = n_elements(input_string_array)
    IF num_entries NE num_redshifts+1 THEN print, 'Crap'

    fnu_array = float(input_string_array[1:num_redshifts])
    fnu_array[where(fnu_array EQ 0.)] = 1.e-30

    hmspec.wave[i] = float(input_string_array[0])
    hmspec.energy[i] = 911.65/hmspec.wave[i]   
    hmspec.fnu[i]  = fnu_array

    readf, unit, input_string
    i = i+1

  ENDWHILE

  close, unit
  free_lun, unit

  mwrfits, /create, hmspec, output_filename
  print, 'Wrote '+output_filename

  RETURN
END
