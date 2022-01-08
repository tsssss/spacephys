;  FILE:
;       d_obj.pro
;
;  CALLING SEQUENCE: d_obj
;
;  PURPOSE:
;       Shows geometric objects (molecules and others.)
;
;  MAJOR TOPICS: Visualizaton and widgets
;
;  CATEGORY:
;       IDL Demo System
;
;  INTERNAL FUNCTIONS and PROCEDURES:
;       fun makespheres         - Create a shere
;       fun makemolecule        - Create a molecule
;       fun moleread            - Read the molecule file
;       pro new_select          - Handle a new selection
;       pro panimate            - Animate (with pattern keyword)
;       pro o3d_animate             - Animate (with pattern keyword)
;       pro toggle_state        - Toggle off and on
;       fun read_noff           - Read the object file
;       pro d_obj_Event       - Event handler
;       pro d_obj_Cleanup    - Cleanup
;       pro d_obj            - Main procedure
;
;  EXTERNAL FUNCTIONS, PROCEDURES, and FILES:
;       pro demo_gettips        - Read the tip file and create widgets
;       pro demo_puttips        - Change tips text
;       object3d.tip
;       knot.nof
;       seashell.nof
;       teapot.nof
;       aspartam.mol
;       caffeine.mol
;       valium.mol
;       pro orb__define.pro  - Create an orb object
;       pro trackball__define.pro  - Create a trackball object
;
;  REFERENCE: IDL Reference Guide, IDL User's Guide
;
;  NAMED STRUCTURES:
;       none.
;
;  COMMON BLOCS:
;       none.
;
;   MODIFICATION HISTORY:
;       9/96,   DD   - Written.
;       10/96,  DAT  - New GUI and combining molecule and objects.
;-

;----------------------------------------------------------------------
;----------------------------------------------------------------------
;
;  Purpose:   toggle off and on state.
;
Function d_objMenuToggleState, wid   ;  IN: widget identifier
    COMPILE_OPT idl2, hidden
    
    WIDGET_CONTROL, wid, GET_VALUE=name
    
    s = STRPOS(name, '(Off)')
    ret = s ne -1                   ;TRUE if new state is on
    if ret then strput, name, '(On )', s $
    else strput, name, '(Off)', strpos(name, '(On )')
    WIDGET_CONTROL, wid, SET_VALUE=name
    RETURN, ret
    end                             ;   of  Toggle_state,
    
    pro d_objMenuSetState, Wid, NAMES=names, Indices = indices, State
    ; Set the vaue of a menu toggle item to State.  Either supply the
    ; Widget ID of the button (wid), or the button's name (wid) along with
    ; the arrays of button names and indices.
    
    COMPILE_OPT idl2, hidden
    
    if n_elements(names) then begin ;Name supplied?
        index = where(strpos(names, Wid) ge 0, count)
        if count le 0 then message, Wid+' not found in menu'
        w = Indices[index[0]]
    endif else w = wid
    
    WIDGET_CONTROL, w, GET_VALUE=value
    On = STRPOS(value, '(On )')
    if (on ge 0) eq State then return ;Already at proper value?
    
    if On ge 0 then strput, value, '(Off)', on $ ;Change value of string
    else strput, value, '(On )', strpos(value, '(Off)')
    
    WIDGET_CONTROL, w, SET_VALUE=value ;update the widget
end


;----------------------------------------------------------------------------
;
;    PURPOSE   Given a uservalue from a menu button created
;              by d_objMenuCreate, return the index of the choice
;              within the category.  Set the selected menu button
;              to insensitive to signify selection, and set all
;              other choices for the category to sensitive.
;
function d_objMapMenuChoice, $
            Eventval, $         ; IN: uservalue from seleted menu button
            MenuItems, $        ; IN: menu item array, as returned by d_objMenuCreate
            MenuButtons         ; IN: button array as returned by d_objMenuCreate

    COMPILE_OPT idl2, hidden
    
    i = STRPOS(eventval, '|', 0)    ;Get the name less the last qualifier
    while (i GE 0) do begin
        j = i
        i = STRPOS(eventval, '|', i+1)
    endwhile
    
    base = STRMID(eventval, 0, j+1) ;  Get the common buttons, includes last | .
    buttons = WHERE(STRPOS(MenuItems, base) EQ 0) ;  buttons that share base name.
    this = (WHERE(eventval EQ MenuItems))[0] ;  Get the Index of the selected item.
    for i=0, N_ELEMENTS(buttons)-1 do begin ;Each button in this category
        index = buttons[i]
        WIDGET_CONTROL, MenuButtons[buttons[i]], $
          SENSITIVE=index NE this
    endfor
    
    RETURN, this - buttons[0]       ;  Return the selected button's index.
end

;----------------------------------------------------------------------------
;
;    PURPOSE  Create a menu from a string descriptor (MenuItems).
;             Return the parsed menu items in MenuItems (overwritten),
;             and the array of corresponding menu buttons in MenuButtons.
;
;    MenuItems = (input/output), on input the menu structure
;                in the form of a string array.  Each button
;                is an element, encoded as follows:
;
;    Character 1 = integer bit flag.  Bit 0 = 1 to denote a
;                  button with children.  Bit 1 = 2 to denote
;                  this is the last child of its parent.
;                  Bit 2 = 4 to show that this button should
;                  initially be insensitive, to denote selection.
;                  Any combination of bits may be set.
;              On RETURN, MenuItems contains the fully
;                  qualified button names.
;
;    Characters 2-end = Menu button text.  Text should NOT
;                       contain the character |, which is used
;                       to delimit menu names.
;
;    Example:
;
;        MenuItems = ['1File', '0Save', '2Quit', $
;       '1Edit', '3Cut', $
;       '3Help']
;
;         Creates a menu with three top level buttons
;         (file, edit and help). File has 2 choices
;         (save and exit), Edit has one choice, and help has none.
;         On RETURN, MenuItems contains the fully qualified
;         menu button names in a string array of the
;         form: ['<Prefix>|File', '<Prefix>|File|Save',
;           '<Prefix>|File|Quit', '<Prefix>|Edit',..., etc. ]
;
pro d_objMenuCreate, $
    MenuItems, $    ; IN/OUT: See below
    MenuButtons, $  ; OUT: Button widget id's of the created menu
    Bar_base, $     ; IN: menu base ID
    Prefix=prefix   ; IN: (opt) Prefix for this menu's button names.
                    ;     If omitted, no prefix

    COMPILE_OPT idl2, hidden
    
    level = 0
    parent = [ bar_base, 0, 0, 0, 0, 0]
    names = STRARR(5)
    lflags = INTARR(5)
    
    MenuButtons = LONARR(N_ELEMENTS(MenuItems))
    
    if (N_ELEMENTS(prefix)) then begin
        names[0] = prefix + '|'
    endif else begin
        names[0] = '|'
    endelse
    
    for i=0, N_ELEMENTS(MenuItems)-1 do begin
        flag = FIX(STRMID(MenuItems[i], 0, 1))
        txt = STRMID(MenuItems[i], 1, 100)
        uv = ''
    
        for j = 0, level do uv = uv + names[j]
        MenuItems[i] = uv + txt     ;  Create the button for fully qualifid names.
        isHelp = 0;txt eq 'Help' or txt eq 'About'
        MenuButtons[i] = WIDGET_BUTTON(parent[level], $
                                       VALUE= txt, UVALUE=uv+txt, $
                                       MENU=flag and 1, HELP=isHelp)
    
        if ((flag AND 4) NE 0) then begin
            WIDGET_CONTROL, MenuButtons[i], SENSITIVE = 0
        endif
    
        if (flag AND 1) then begin
            level = level + 1
            parent[level] = MenuButtons[i]
            names[level] = txt + '|'
            lflags[level] = (flag and 2) NE 0
        endif else if ((flag AND 2) NE 0) then begin
            while lflags[level] do level = level-1 ;  Pops the previous levels.
            level = level - 1
        endif
    endfor
end


;----------------------------------------------------------------------
;
;  Purpose: This reads a modified .off format
;           file (header removed, one ptype...)
;
function d_objReadNoff, $
    file, $            ; IN: filename
    xr, $              ; OUT: x radius
    yr, $              ; OUT: y radius
    zr                 ; OUT: z radius

    COMPILE_OPT idl2, hidden
    
    t0 = systime(1)
    s = ' '
    npsize = 1
    
    RESTORE, file
    
    xr = [min(x, max = xx), xx]     ;Get ranges
    yr = [min(y, max = xx), xx]
    zr = [min(z, max = xx), xx]
    
    sc = [xr[1]-xr[0], yr[1]-yr[0], zr[1]-zr[0]] ;Ranges
    xr[0] = (xr[1] + xr[0])/2.0     ;Midpoint
    yr[0] = (yr[1] + yr[0])/2.0
    zr[0] = (zr[1] + zr[0])/2.0
    s = max(sc)                     ;Largest range...
    
    x = (x - xr[0]) / s
    y = (y - yr[0]) / s
    z = (z - zr[0]) / s
    
    xr = [-0.7, 0.7]                ;Fudge the ranges
    yr = xr
    zr = xr
    s = OBJ_NEW("IDLgrPolygon", TRANSPOSE([[x],[y],[z]]), $
                SHADING=1, $
                POLY=mesh, COLOR=[200,200,200])
    
    ; print, 'D_OBJREADNOFF:', systime(1)-t0
    RETURN, s
    
end                             ;   of d_objReadNoff


;----------------------------------------------------------------------
;
;  Purpose:  Read the molecule data input file.
;
function d_objMolRead, $
    filename, $          ; IN: filename
    atom_xyz, $
    atoms

    COMPILE_OPT idl2, hidden
    
    n_atoms = 0L
    atom_xyz = FLTARR(3,1000, /nozero)
    x = (y = (z = 0.))
    atoms = INTARR(1000, /nozero)
    
    OPENR, lun, filename, /GET_LUN, error = i
    if i lt 0 then message, filename + ' not found'
    s = " "
    
    while eof(lun) eq 0 do begin
        READF, lun, s
        i = strpos(s, ":")
        if i NE -1 then begin
            READS, strmid(s,i+1,strlen(s)-i+1),x,y,z
            case strmid(s,0,2) of
                "C:" : a_type = 0
                "H:" : a_type = 1
                "N:" : a_type = 2
                "S:" : a_type = 3
                "X:" : a_type = 4
                "Ce" : a_type = 5
                "BR" : a_type = 6
                "O:" : a_type = 7
                else : a_type = -1
            endcase
    
            if (a_type NE -1) then begin
                atom_xyz[0,n_atoms] = x
                atom_xyz[1,n_atoms] = y
                atom_xyz[2,n_atoms] = z
                atoms[n_atoms] = a_type
                n_atoms = n_atoms + 1
                if n_atoms ge n_elements(atoms) then $
                  message, 'Too many atoms in molecule'
            endif
        endif
    endwhile
    
    endoffile:
    
    FREE_LUN,lun
    
    if (n_atoms GE 0) then begin
        atoms = atoms[0:n_atoms-1]
        atom_xyz = atom_xyz[*,0:n_atoms-1]
    endif
    
    RETURN, n_atoms
end   ;   of  d_objMolRead




function d_objMakeSpheres, DENSITY=dens ;Make spheres for the 8 types of atoms
    COMPILE_OPT idl2, hidden
    
    o = objarr(8)
    if n_elements(dens) eq 0 then dens = .9
    
    o[0] = OBJ_NEW('orb', COLOR=[128,128,128], DENSITY=dens, RADIUS=1.5) ; C
    o[1] = OBJ_NEW('orb', COLOR=[255,255,255], DENSITY=dens, RADIUS=0.9) ; H
    o[2] = OBJ_NEW('orb', COLOR=[0,0,255],     DENSITY=dens, RADIUS=1.5) ; N
    o[3] = OBJ_NEW('orb', COLOR=[230,230,0],   DENSITY=dens, RADIUS=1.5) ; S
    o[4] = OBJ_NEW('orb', COLOR=[25,25,0],     DENSITY=dens, RADIUS=1.5) ; X
    o[5] = OBJ_NEW('orb', COLOR=[204,0,153],   DENSITY=dens, RADIUS=4.7) ; Ce
    o[6] = OBJ_NEW('orb', COLOR=[204,25,230],  DENSITY=dens, RADIUS=1.3) ; BR
    o[7] = OBJ_NEW('orb', COLOR=[255,0,0],     DENSITY=dens, RADIUS=1.46) ; O
    RETURN, o
end



Pro d_objAddLight, Model             ;Set up the lighting
    COMPILE_OPT idl2, hidden
    
    Model->Add, OBJ_NEW('IDLgrLight', $ ;Directional
                        LOCATION=[2,2,5], TYPE=2, COLOR=[255,255,255], $
                        INTENSITY=0.6 )
    
    Model->Add, OBJ_NEW('IDLgrLight', $ ;Directional
                        LOCATION=[2,-2,-5], TYPE=2, COLOR=[255,255,0], $
                        INTENSITY=0.25 )
    
    Model->add, OBJ_NEW('IDLgrLight', TYPE=0, $ ;Ambient light
                        INTENSITY=0.375, COLOR=[255,255,255])
end


Function d_objMakeMolecule, filename, state
    COMPILE_OPT idl2, hidden
    
    t0 = systime(1)
    if state.spheres[0] eq obj_new() then begin ;Init spheres?
        state.spheres = d_objMakeSpheres(DENSITY=0.8) ;The spheres for the atoms
    endif
    
    n = d_objMolRead(filename, atom_xyz, atom_type) ;Get atoms
    sc = 0.3 * 2.0 / (max(atom_xyz[0,*], min = mn) - mn) ;X extent
    
    state.oModelTop = OBJ_NEW('IDLgrModel')       ;Top model
    state.oModelSurface = OBJ_NEW('IDLgrModel')
    state.oModelTop->add, state.oModelSurface
    ;state.oModelSurface->add, state.oSurface
    
    state.oModelOffset = OBJ_NEW('IDLgrModel')
    state.oModelOffset->translate, 0, 0, 0.005   ;Offset Z to make visible
    state.oModelEdges = OBJ_NEW('IDLgrModel')
    state.oModelOffset->add, state.oModelEdges
    
    
    ; destroy the old objects prior to saving the new ones
    obj_destroy, state.oSurface
    state.oSurface = OBJ_NEW('IDL_CONTAINER')
    names = ['Carbon','Hydrogen','Nitrogen','Silicon', 'Xenon', $
             'Cesium', 'Boron', 'Oxygen']
    
    for i=0,n-1 do begin            ;For each atom in molecule
        s = state.spheres[atom_type[i]]   ;A sphere for the atom
        s->GetProperty, POBJ=sh
        sh->GetProperty, COLOR=col, POLY=pmesh
        p = OBJ_NEW('IDLgrPolygon', SHARE_DATA=sh, POLY=pmesh, COLOR=col)
        p->SetProperty, SHADING=1
        oModelAtom = OBJ_NEW('IDLgrModel', /SELECT_TARGET)
        oModelAtom->SetProperty, UVALUE=names[atom_type[i] > 0 < (n_elements(names)-1)]
        oModelAtom->Translate, atom_xyz[0,i], atom_xyz[1,i], atom_xyz[2,i]
        oModelAtom->Add, p
        state.oSurface->add, p
        state.oModelSurface->Add, oModelAtom
        ;state.ModelTop->Add, state.oModelSurface
                                    ; Make the mesh object
        oModelAtomEdge = OBJ_NEW('IDLgrModel')
        p1 = OBJ_NEW('IDLgrPolyline', SHARE_DATA=sh, POLYLINE=pmesh, COLOR=[0,0,0])
    ;    p1->SetProperty, SHADING=1
        ;state.oModelEdges = OBJ_NEW('IDLgrModel')
        state.oModelEdges->add, oModelAtomEdge
        oModelAtomEdge->add, p1
        ; Don't add offset here.  Do it below to parent model
        oModelAtomEdge->translate, atom_xyz[0,i], $
            atom_xyz[1,i], atom_xyz[2,i]
        ;state.oModelEdges->add, state.oModelEdges
    endfor
    
    state.oModelOffset->translate, 0, 0, 0.005
    state.oModelSurface->Scale, sc, sc, sc
    state.oModelEdges->Scale, sc, sc, sc
    
    d_objAddLight, state.oModelTop
    
    if state.vert_coloring then begin
        state.vert_coloring = 0
        d_objMenuSetState, '|Options|Vertex Coloring', 0, $
          Names=state.menuitems, Indices=state.menubuttons
    endif
    
    state.oModelTop->add, state.oModelOffset
    
    demo_puttips, state, ['molec', 'right'], [10,11], /LABEL
    
    ; print, 'd_objMakeMolecule:', systime(1) - t0
    return, state.oModelTop
end



Function d_objMakeShape, filename, state
    COMPILE_OPT idl2, hidden
    
    state.oModelTop = OBJ_NEW('IDLgrModel')       ;Top model
    
    state.oSurface = d_objReadNoff(filename, xr, yr, zr)
    
    state.oModelSurface = OBJ_NEW('IDLgrModel')
    state.oModelTop->add, state.oModelSurface
    state.oModelSurface->add, state.oSurface
    
    state.oModelOffset = OBJ_NEW('IDLgrModel')
    state.oModelOffset->translate, 0, 0, 0.005   ;Offset Z to make visible
    state.oModelEdges = OBJ_NEW('IDLgrModel')
    state.oModelOffset->add, state.oModelEdges
    state.oSurface->GetProperty, POLY=pmesh
    p1 = OBJ_NEW('IDLgrPolyline', SHARE_DATA=state.oSurface, POLY=pmesh, $
       COLOR=[0,0,0])
    state.oModelEdges->add, p1                  ;Add the edging data
    
;    if state.vert_coloring eq 0 then begin ;Default = vertex coloring ON
;        state.vert_coloring = 1
;        d_objMenuSetState, '|Options|Vertex Coloring', 1, $
;          NAMES=state.menuitems, INDICES=state.menubuttons
;    endif
    
    state.oModelTop->add, state.oModelOffset
    
    d_objAddLight, state.oModelTop
    
;    demo_puttips, state, ['inter', 'mouse'], [10,11], /LABEL
    
    return, state.oModelTop
end




pro d_objSetObjectAttributes, state ;Set object's attrib to current settings
    COMPILE_OPT idl2, hidden
    
    if state.oModelEdges ne obj_new() then $
      d_objSetProp, state.oModelEdges, HIDE=1-state.edging
    
    if state.vert_coloring then vc = state.vc else vc = 0
    if state.bottomColor then bot = 0 else bot = [64, 192, 128]
    
    d_objSetProp, state.oSurface, STYLE=([0,1,2,0,1])[state.style], $
      HIDDEN_LINE=([0,0,0,1,1])[state.style], BOTTOM=bot, REJECT = state.reject, $
      THICK=state.thick, VERT_COLORS=vc, SHADING=state.shading
end


pro d_objLoadItem, index, State      ;Load the new item...
    COMPILE_OPT idl2, hidden
    
    filename = demo_filepath(state.items[2,index], $ ;Input file name
                        SUBDIR=['examples','demo', 'demodata'])
    
    
    state.oModelTop = call_function(state.items[1,index], filename, state) ;Load it
    
    
    state.oModelSurface->GetProperty, TRANSFORM=tmg0
    ;state.omodel = o
    
    state.oView->add, state.oModelTop
    
    state.tmg0 = tmg0
    ;o1 = state.oSurface
    
    ;if obj_class(o1) eq 'IDL_CONTAINER' then o1 = o1->get(/all)
    ;for i=0, n_elements(o1)-1 do begin
    ;    o1(i)->GetProperty, POLY=pmesh
    ;    p = OBJ_NEW('IDLgrPolyline', COLOR=[0,0,0], $
    ;                SHARE_DATA = o1(i), POLYLINE=pmesh, HIDE=1)
    ;    state.oModelEdges->add, p
    ;    state.omodel->add, p
    ;endfor
    
    ; if state.oModelEdges ne obj_new() then
    ;state.omodel->add, state.oModelEdges
    ;state.oModelTop->add, state.oModelOffset
    state.cur_sel = obj_new()       ;Nothing's selected
    d_objSetObjectAttributes, state    ; Now set the attributes to current values
end



PRO d_objSetProp, o, _extra = e ;Set a property on either a container that
    ; contains a number of objects, or a single object.
    
    COMPILE_OPT idl2, hidden
    
    if obj_class(o) eq 'IDL_CONTAINER' then arr = o->get(/all) $
    else arr = o
    for i=0, n_elements(arr)-1 do arr[i]->setproperty, _EXTRA=e
end



pro d_objNewSelect, state, target
    COMPILE_OPT idl2, hidden
    
    ostyle = ([0,1,2,0,1])[state.style] ;Original style
    if (state.cur_sel NE obj_new()) then begin ;Remove old selection
        state.cur_sel->SetProperty,STYLE=ostyle
        state.cur_sel = obj_new()
    endif
    
    if (N_ELEMENTS(target) NE 0) then begin ;Add new selection
        if ostyle eq 2 then nstyle = 1 else nstyle = 2 ;New style
        state.cur_sel = target
        state.cur_sel->SetProperty, STYLE=nstyle
    endif
end                             ;     of d_objNewSelect


pro d_objDrawEvent, ev        ;Handle events for the draw window
    COMPILE_OPT idl2, hidden
    
    WIDGET_CONTROL, ev.top, GET_UVALUE=state, /NO_COPY
    
    if (ev.type EQ 4) then begin ;  Expose.
        WIDGET_CONTROL, ev.top, /HOURGLASS ;Redraw entire window
        state.oWindow->draw, state.oView
    
    endif else begin                ;Not expose
        bHaveTransform = state.oTrack->Update(ev, TRANSFORM=qmat) ;trackball update
        if (bHaveTransform NE 0) then begin
            state.oModelSurface->GetProperty, TRANSFORM=t
            state.oModelSurface->SetProperty,TRANSFORM= t # qmat
            state.oModelEdges->SetProperty,TRANSFORM= t # qmat
        endif
    
        if (ev.type EQ 0) then begin ;  Button press.
            if ev.press gt 1 then begin ;Press with right or middle button
                picked = state.oWindow->select(state.oView, [ev.x, ev.y])
                if (size(picked))[0] eq 0 then goto, done ;Hit anything?
                if obj_class(picked[0]) ne "IDLGRMODEL" then goto, done
                picked[0]->GetProperty, UVALUE=uval, TRANSFORM=tm
                str = string(uval, tm[3,0], tm[3,1], tm[3,2], $
                             format='(A, " at [", 3F6.2, "]")')
                d_objNewSelect, state, picked[0]->get()
                demo_puttips, state, [str, ''], [10,11] ;Label the molecule
                goto, draw_it
            endif else begin        ;Not rt or middle button
                state.btndown = 1
                state.oWindow->SetProperty, QUALITY=state.dragq
                WIDGET_CONTROL, state.wDraw, /DRAW_MOTION
            endelse
        endif else if ((ev.type eq 2) and (state.btndown eq 1)) then begin ;MOTION.
            if (bHaveTransform) then state.oWindow->Draw, state.oView
        endif else if (ev.type eq 1) and state.btndown then begin ;Release
            state.btndown = 0
            state.oWindow->SetProperty, QUALITY=2
            draw_it:
            WIDGET_CONTROL, ev.top, /HOURGLASS
            t0 = systime(1)
            state.oWindow->Draw, state.oView
            ;demo_puttips, state, $
            ;  'Time =' + STRING(systime(1)-t0, FORMAT='(F6.2)')+ ' seconds', 12
            WIDGET_CONTROL, state.wDraw, DRAW_MOTION=0
        endif
    endelse
    
    done: WIDGET_CONTROL, ev.top, SET_UVALUE=state, /NO_COPY
    
end




pro d_objEvent, ev             ;Main event handler
    COMPILE_OPT idl2, hidden
    
    if (TAG_NAMES(ev, /STRUCTURE_NAME) EQ $ ;  Quit from the close box.
        'WIDGET_KILL_REQUEST') then begin
        WIDGET_CONTROL, ev.top, /DESTROY
        RETURN
    endif
    
    WIDGET_CONTROL, ev.id, GET_UVALUE=uval
    WIDGET_CONTROL, ev.top, GET_UVALUE=state, /NO_COPY
    
    
    if uval eq '|File|Quit' then  begin
        WIDGET_CONTROL, ev.top, /DESTROY, SET_UVALUE=state, /NO_COPY
        return
    
    endif else if uval eq 'SCALING' then begin
        state.scaling = ev.value / 100. ;New scale factor
        goto, new_scale
    
    endif else if uval eq 'RESET' then begin ;
        if ev.value eq 'Defaults' then begin ;Reset the defaults?
            state.shading = 1
            state.thick = 1
            state.edging = 0
            state.reject = 0
            state.bottomColor = 1
            state.style = 2
            d_objMenuSetState, '|Edging', state.edging, $
              Names=state.menuitems, Indices=state.menubuttons
            d_objMenuSetState, '|Bottom Color', state.bottomColor, $
              Names=state.menuitems, Indices=state.menubuttons
            d_objMenuSetState, '|Backface', state.reject, $
              Names=state.menuitems, Indices=state.menubuttons
        endif
    
    reset_transform:
        state.scaling = 1.0
        WIDGET_CONTROL, state.wScalingSlider, SET_VALUE=100 ;And scaling
        state.oModelSurface->SetProperty, TRANSFORM=state.tmg0 ;And fall thru
        state.oModelEdges->SetProperty, TRANSFORM=state.tmg0 ;And fall thru
    new_scale:
        s = state.scaling
        ;state.oModelSurface->SetProperty, $
        ;  TRANSFORM = [[s,0,0,0],[0,s,0,0], [0,0,s,0], [0,0,0,1]] ;New transform
        ;state.oModelEdges->SetProperty, $
        ;  TRANSFORM = [[s,0,0,0],[0,s,0,0], [0,0,s,0], [0,0,0,1]] ;New transform
        state.oModelTop->SetProperty, $
          TRANSFORM = [[s,0,0,0],[0,s,0,0], [0,0,s,0], [0,0,0,1]] ;New transform
    
    redraw:
        d_objSetObjectAttributes, state ;Update attributes
        widget_control, ev.top, /HOURGLASS
        t0 = systime(1)
        state.oWindow->Draw, state.oView
        ;demo_puttips, state, $
        ;  'Time =' + STRING(systime(1)-t0, FORMAT='(F6.2)')+ ' seconds', 12
    
    endif else if strpos(uval, '|Shading') ge 0 then begin
        state.shading = d_objMapMenuChoice(uval, state.MenuItems, $
                           state.MenuButtons)
        goto, redraw
    
    endif else if strpos(uval, '|Drag ') ge 0 then begin
        state.dragq = d_objMapMenuChoice(uval, state.MenuItems, state.MenuButtons)
    endif else if strpos(uval, '|Vertex') ge 0 then begin
        state.vert_coloring = d_objMenuToggleState(ev.id)
        goto, redraw
    
    endif else if strpos(uval, '|Line ') ge 0 then begin
        state.thick = ([1,3,5,7])[d_objMapMenuChoice(uval, state.MenuItems, $
                                                  state.MenuButtons)]
        goto, redraw
    
    endif else if strpos(uval, '|Edging') ge 0 then begin
        state.edging = d_objMenuToggleState(ev.id)
        goto, redraw
    
    endif else if strpos(uval, '|Backface') ge 0 then begin
        state.reject = d_objMenuToggleState(ev.id) ;New backface culling
        goto, redraw
    
    endif else if strpos(uval, '|Bottom Color') ge 0 then begin
        state.bottomColor = d_objMenuToggleState(ev.id)    ;New backface culling
        goto, redraw
    
    endif else if strpos(uval, '|Style') ge 0 then begin
        state.style = d_objMapMenuChoice(uval, state.MenuItems, state.MenuButtons)
        goto, redraw
    
    endif else if uval eq 'OBJ' then begin ;Load a new object
        prev = state.oModelTop
        state.oView->remove, prev
        obj_destroy, state.oModelTop
        d_objLoadItem, ev.value, state
        obj_destroy, prev
        goto, reset_transform
    endif else begin
        print, "Dunno: ", uval
    endelse
    
    
    done : WIDGET_CONTROL, ev.top, SET_UVALUE=state, /NO_COPY
end


pro d_objCleanup, wBase
    COMPILE_OPT idl2, hidden
    
    WIDGET_CONTROL, wBase, GET_UVALUE=state, /NO_COPY
    
    OBJ_DESTROY, state.oView        ;  Destroy the top objects
    OBJ_DESTROY, state.oTrack
    OBJ_DESTROY, state.oModelTop
    
    OBJ_DESTROY, state.spheres
    OBJ_DESTROY, state.oSurface
    OBJ_DESTROY, state.OModelEdges
    
    TVLCT, state.colorTable         ;  Restore the color table.
    
    if widget_info(state.groupBase, /valid) then $
      widget_control, state.groupBase, /map
    
end                             ;  of d_objCleanup





PRO d_obj, $
    PATTERN=pattern, $         ; IN: (opt) line pattern
    GROUP=group, $             ; IN: (opt) group identifier
    RECORD_TO_FILENAME=record_to_filename, $
    APPTLB = appTLB            ; OUT: (opt) TLB of this application
    
    COMPILE_OPT idl2
    
    t0 = systime(1)
    
    ; Check the validity of the group identifier
    ;
    ngroup = N_ELEMENTS(group)
    if (ngroup NE 0) then begin
        check = WIDGET_INFO(group, /VALID_ID)
        if (check NE 1) then begin
            print,'Error, the group identifier is not valid'
            print, 'Return to the main application'
            RETURN
        endif
        groupBase = group
    endif else groupBase = 0L
    
    
    
    Device, GET_SCREEN_SIZE=screenSize ;Size the drawable : to the screen size
    xdim = floor(screenSize[0] * 0.6)
    ydim = floor(xdim * 0.8)
    
    TVLCT, savedR, savedG, savedB, /GET ;  save current color table
    colorTable = [[savedR],[savedG],[savedB]]
    
    Items = [ ['Seashell', 'd_objMakeShape', 'seashell.dat'], $
        ['Knot', 'd_objMakeShape', 'knot.dat'], $
        ['Teapot', 'd_objMakeShape', 'teapot.dat'], $
        ['Valium Molecule', 'd_objMakeMolecule', 'valium.mol'], $
        ['Aspartame Molecule', 'd_objMakeMolecule', 'aspartam.mol'], $
        ['Caffeine Molecule', 'd_objMakeMolecule', 'caffeine.mol']]
    nObjects = N_ELEMENTS(Items)/3
    
    if (N_ELEMENTS(group) EQ 0) then begin
       wBase = WIDGET_BASE(/COLUMN, $
            XPAD=0, YPAD=0, $
            /TLB_KILL_REQUEST_EVENTS, $
            TLB_FRAME_ATTR=1, MBAR=barBase, TITLE="Three Dimensional Geometry")
    endif else begin
        wBase = WIDGET_BASE(/COLUMN, $
            XPAD=0, YPAD=0, $
            /TLB_KILL_REQUEST_EVENTS, $
            GROUP_LEADER=group, $
            TLB_FRAME_ATTR=1, MBAR=barBase, TITLE="Three Dimensional Geometry")
    endelse

    wDraw = widget_draw(wBase, GRAPHICS_LEVEL=2, $
        XSIZE=xdim, YSIZE=ydim, /BUTTON_EVENTS, $
        UVALUE='DRAW', RETAIN=0, /EXPOSE_EVENTS, $
        EVENT_PRO = 'd_objDrawEvent')
    
    wStatusBase = WIDGET_BASE(wBase, /ROW) ;Tips base
    
    WIDGET_CONTROL, wBase, /REALIZE
    
    appTLB = wBase    ;  Return the top level base in the appTLB keyword.
    
    WIDGET_CONTROL, wDraw, GET_VALUE=oWindow ; Grab the window id of the drawable.
    
;    sText = demo_getTips(demo_filepath("object3d.tip", $
;        SUBDIR=['examples','demo','demotext']), $
;        wBase, wStatusBase) ;Get the tips
    
    aspect = float(xdim)/float(ydim) ; viewplane rect based on aspect ratio.
    myview = [-0.5,-0.5,1,1]
    if (aspect > 1) then begin
        myview[0] = myview[0] - ((aspect-1.0)*myview[2])/2.0
        myview[2] = myview[2] * aspect
    endif else begin
        myview[1] = myview[1] - (((1.0/aspect)-1.0)*myview[3])/2.0
        myview[3] = myview[3] * aspect
    endelse
    
    oView = OBJ_NEW('idlgrview', PROJECTION=2, EYE=3, ZCLIP=[1.5,-1.5],$
        VIEWPLANE_RECT=myview, COLOR=[90,90,90])
    
    oTrack = OBJ_NEW('Trackball', [xdim/2.0, ydim/2.0], xdim/2.0)
    
    
    
    state = $
      { wBase: wBase, $             ; Main base
        wDraw: wDraw, $             ; Widget draw ID
        groupBase: groupBase, $     ; Base of Group Leader
        spheres : objarr(8), $      ;spheres for molecules
        currentItem : 0, $          ;Index of current item
        items : items, $            ;array of object items... (3,nitems)
        oView : oView, $            ;the view
        oModelTop: obj_new(), $
        oModelSurface: obj_new(), $
        oModelEdges: obj_new(), $
        oModelOffset: obj_new(), $
        oTrack: oTrack, $
        oWindow : oWindow, $
        oSurface : obj_new(), $     ;Visible surface
        tmg0: fltarr(4,4), $        ;Initial transform
        dragq: 1, $                 ;Drag quality
        btndown: 0, $               ;Button down flag
        colorTable: colorTable, $
        cur_sel: obj_new(), $
        scaling: 1.0, $
        shading: 1, $
        vert_coloring: 0, $
        thick: 1, $
        edging: 0, $
        reject: 0, $
        bottomColor: 1, $
        style: 2 $
      }
    
    d_objLoadItem, 0, state              ;Load the first item
    WIDGET_CONTROL, wBase, SET_UVALUE=state, /NO_COPY
    
    ;print, systime(1) - t0, ' seconds'
    XMANAGER, 'd_obj', wBase, /NO_BLOCK, $
       EVENT_HANDLER='d_objEvent', $
       CLEANUP='d_objCleanup'
end   ;   of   d_obj
