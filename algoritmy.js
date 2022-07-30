/*   **************************************************************************************************
tato funkce BABYLON.Mesh.prototype.increaseFacets je převzatá od tretí strany
je dostupná na adrese https://doc.babylonjs.com/toolsAndResources/utilities/Increasing_Facets
slouží k zvýšení počtu FaceT a tak i ke tvýšení přesnosti algoritmů, které FaceT využívají
*/

    BABYLON.Mesh.prototype.increaseFacets = function(pps) { 
    var _gaps = pps+1;
    var _n = _gaps + 1;
    var _fvs =[];
    for(var _i=0; _i<_n; _i++) {
        _fvs[_i] = [];
    }    
    var _A,_B;
    var _d ={x:0,y:0,z:0};
    var _u ={x:0,y:0};
    var _indices = [];
    var _vertexIndex = [];
    var _side = [];
	var _l; 
    var _uvs = this.getVerticesData(BABYLON.VertexBuffer.UVKind);
    var _meshIndices = this.getIndices();
    var _positions = this.getVerticesData(BABYLON.VertexBuffer.PositionKind);    
    var _normals =[];    

    for(var _i = 0; _i<_meshIndices.length; _i+=3) {
        _vertexIndex[0] = _meshIndices[_i];
        _vertexIndex[1] = _meshIndices[_i + 1];
        _vertexIndex[2] = _meshIndices[_i + 2];        
        for(var _j = 0; _j<3; _j++) {
            _A = _vertexIndex[_j];
            _B = _vertexIndex[(_j+1)%3];        
            if(_side[_A] === undefined  && _side[_B] ===  undefined) {            
                _side[_A] = [];
                _side[_B] = [];            
            }
            else {
                if(_side[_A] === undefined) {                    
                    _side[_A] = [];
                }
                if(_side[_B] === undefined) {                    
                    _side[_B] = [];                                
                }
            }
            if(_side[_A][_B]  === undefined  && _side[_B][_A] === undefined) {            
                _side[_A][_B] = [];
                _d.x = (_positions[3 * _B] - _positions[3 * _A])/_gaps;
                _d.y = (_positions[3 * _B + 1] - _positions[3 * _A + 1])/_gaps;
                _d.z = (_positions[3 * _B + 2] - _positions[3 * _A + 2])/_gaps;
                _u.x = (_uvs[2*_B] - _uvs[2*_A])/_gaps;
                _u.y = (_uvs[2*_B + 1] - _uvs[2*_A + 1])/_gaps;
                _side[_A][_B].push(_A);                
                for(var _k=1; _k<_gaps; _k++) {                
                    _side[_A][_B].push(_positions.length/3);                
                    _positions.push(_positions[3 * _A] + _k*_d.x, _positions[3 * _A + 1] + _k*_d.y, _positions[3 * _A + 2] + _k*_d.z);
                    _uvs.push(_uvs[2*_A] + _k*_u.x, _uvs[2*_A + 1] + _k*_u.y);
                }                
                _side[_A][_B].push(_B);
                _side[_B][_A]=[];
                _l = _side[_A][_B].length;
                for(var _a=0; _a<_l; _a++) {
                    _side[_B][_A][_a] = _side[_A][_B][_l-1-_a];
                }
            }
            else {
                if(_side[_A][_B] === undefined) {            
                    _side[_A][_B]=[];
                    _l = _side[_B][_A].length;
                    for(var _a=0; _a<_l; _a++) {
                        _side[_A][_B][_a] = _side[_B][_A][_l-1-_a];
                    }
                }
                if(_side[_B][_A] === undefined) {            
                    _side[_B][_A]=[];                
                    _l = _side[_A][_B].length;
                    for(var _a=0; _a<_l; _a++) {
                        _side[_B][_A][_a] = _side[_A][_B][_l-1-_a];
                    }
                }
            }                    
        }    
        _fvs[0][0] = _meshIndices[_i];
        _fvs[1][0] = _side[_meshIndices[_i]][_meshIndices[_i + 1]][1];
        _fvs[1][1] = _side[_meshIndices[_i]][_meshIndices[_i + 2]][1];        
        for(var _k = 2; _k<_gaps; _k++) {
            _fvs[_k][0] = _side[_meshIndices[_i]][_meshIndices[_i + 1]][_k];
            _fvs[_k][_k] = _side[_meshIndices[_i]][_meshIndices[_i + 2]][_k];        
            _d.x = (_positions[3 * _fvs[_k][_k]] - _positions[3 * _fvs[_k][0]])/_k;
            _d.y = (_positions[3 * _fvs[_k][_k] + 1] - _positions[3 * _fvs[_k][0] + 1])/_k;
            _d.z = (_positions[3 * _fvs[_k][_k] + 2] - _positions[3 * _fvs[_k][0] + 2])/_k;
            _u.x = (_uvs[2*_fvs[_k][_k]] - _uvs[2*_fvs[_k][0]])/_k;
            _u.y = (_uvs[2*_fvs[_k][_k] + 1] - _uvs[2*_fvs[_k][0] + 1])/_k;
            for(var _j = 1; _j<_k; _j++) {                
                _fvs[_k][_j] = _positions.length/3;                
                _positions.push(_positions[3 * _fvs[_k][0]] + _j*_d.x, _positions[3 * _fvs[_k][0] + 1] + _j*_d.y, _positions[3 * _fvs[_k][0] + 2] + _j*_d.z);
                _uvs.push(_uvs[2*_fvs[_k][0]] + _j*_u.x, _uvs[2*_fvs[_k][0] + 1] + _j*_u.y);
            }        
        }
        _fvs[_gaps] = _side[_meshIndices[_i + 1]][_meshIndices[_i + 2]];

        _indices.push(_fvs[0][0],_fvs[1][0],_fvs[1][1]);
        for(var _k = 1; _k<_gaps; _k++) {
            for(var _j = 0; _j<_k; _j++) {            
                _indices.push(_fvs[_k][_j],_fvs[_k+1][_j],_fvs[_k+1][_j+1]);
                _indices.push(_fvs[_k][_j],_fvs[_k+1][_j+1],_fvs[_k][_j+1]);
            }        
            _indices.push(_fvs[_k][_j],_fvs[_k+1][_j],_fvs[_k+1][_j+1]);
        }

    }                            

    var vertexData = new BABYLON.VertexData();
    vertexData.positions = _positions;
    vertexData.indices = _indices;
    vertexData.uvs = _uvs;

    BABYLON.VertexData.ComputeNormals(_positions, _indices, _normals);
    vertexData.normals = _normals;
	
    vertexData.applyToMesh(this);

    }  // konec prevzate funkce
    
    /*
    funkce, která zjistí zda jsou objekty dostatečně blízko sebe aby vůbec mělo smysl objekty
    zkoumat na kolizi pomalým a presným algortimem
    využíva poznatky o AABB
    autor: Tomáš Pribyl
    studijní materiály: [1][5][6][7]
    */
    function optimalizacni_funkce(obj1, obj2)
    {
        var posx = obj1.position.x - obj2.position.x;                // zjisteni vzdalenosti objektu
            var posy = obj1.position.y - obj2.position.y;
            var posz = obj1.position.z - obj2.position.z;
            
            if(posx<0){posx = posx*-1}
            if(posy<0){posy = posy*-1}
            if(posz<0){posz = posz*-1}
            
            var rozmer_herox = obj1.getBoundingInfo().maximum["x"] ;     // zjisteni rozmeru objektu
            var rozmer_heroy = obj1.getBoundingInfo().maximum["y"] ;
            var rozmer_heroz = obj1.getBoundingInfo().maximum["z"] ;
            
            var rozmer_groundx = obj2.getBoundingInfo().maximum["x"];
            var rozmer_groundy = obj2.getBoundingInfo().maximum["y"];
            var rozmer_groundz = obj2.getBoundingInfo().maximum["z"];
            
            //xxxxxxxxxxxxxxxxxxxxxxxxxxx                                       // objekty se nyni zrotuji a vypocita se jejich novy rozmer
            var euler_ground = obj2.rotationQuaternion.toEulerAngles();        // rotace objektu v ose x
            
            var pripona_ground = Math.sqrt(rozmer_groundx*rozmer_groundx  +  rozmer_groundy*rozmer_groundy);
            
            var uhelx_ground = Math.asin(rozmer_groundx/pripona_ground);
            
            rozmer_groundx = pripona_ground * Math.sin(euler_ground.z + uhelx_ground);
            rozmer_groundy = pripona_ground * Math.cos(euler_ground.z + uhelx_ground);
            
            if(rozmer_groundx<0){rozmer_groundx=rozmer_groundx*-1}
            if(rozmer_groundy<0){rozmer_groundy=rozmer_groundy*-1}
            
            //xxxxxxxxxxxxxxxxxxxxxxx
            
            pripona_ground = Math.sqrt(rozmer_groundx*rozmer_groundx  +  rozmer_groundz*rozmer_groundz);      // rotace objektu v ose z
            
            var uhelz_ground = Math.asin(rozmer_groundz/pripona_ground);
            
            rozmer_groundz = pripona_ground * Math.sin(euler_ground.y + uhelz_ground);
            rozmer_groundx = pripona_ground * Math.cos(euler_ground.y + uhelz_ground);
            
            if(rozmer_groundx<0){rozmer_groundx=rozmer_groundx*-1}
            if(rozmer_groundz<0){rozmer_groundz=rozmer_groundz*-1}
            
            //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            
            pripona_ground = Math.sqrt(rozmer_groundz*rozmer_groundz  +  rozmer_groundy*rozmer_groundy);   // rotace objektu v ose y
            
            var uhely_ground = Math.asin(rozmer_groundy/pripona_ground);
            
            rozmer_groundy = pripona_ground * Math.sin(euler_ground.x + uhely_ground);
            rozmer_groundz = pripona_ground * Math.cos(euler_ground.x + uhely_ground);
            
            if(rozmer_groundy<0){rozmer_groundy=rozmer_groundy*-1}
            if(rozmer_groundz<0){rozmer_groundz=rozmer_groundz*-1}
            
            //yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
            var euler_hero = obj1.rotationQuaternion.toEulerAngles()             // rotace objektu v ose x
            var pripona_hero = Math.sqrt(rozmer_herox*rozmer_herox  +  rozmer_heroy*rozmer_heroy);
            var uhelx_hero = Math.asin(rozmer_herox/pripona_hero);
            rozmer_herox = pripona_hero * Math.sin(euler_hero.z + uhelx_hero);
            rozmer_heroy = pripona_hero * Math.cos(euler_hero.z + uhelx_hero);
            if(rozmer_herox<0){rozmer_herox=rozmer_herox*-1}
            if(rozmer_heroy<0){rozmer_heroy=rozmer_heroy*-1}
            //yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
            pripona_hero = Math.sqrt(rozmer_herox*rozmer_herox  +  rozmer_heroz*rozmer_heroz);       // rotace objektu v ose z
            var uhelz_hero = Math.asin(rozmer_heroz/pripona_hero);
            rozmer_heroz = pripona_hero * Math.sin(euler_hero.y + uhelz_hero);
            rozmer_herox = pripona_hero * Math.cos(euler_hero.y + uhelz_hero);
            if(rozmer_herox<0){rozmer_herox=rozmer_herox*-1}
            if(rozmer_heroz<0){rozmer_heroz=rozmer_heroz*-1}
            //yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
            pripona_hero = Math.sqrt(rozmer_heroz*rozmer_heroz  +  rozmer_heroy*rozmer_heroy);       // rotace objektu v ose y
            var uhely_hero = Math.asin(rozmer_heroy/pripona_hero);
            rozmer_heroy = pripona_hero * Math.sin(euler_hero.x + uhely_hero);
            rozmer_heroz = pripona_hero * Math.cos(euler_hero.x + uhely_hero);
            if(rozmer_heroy<0){rozmer_heroy=rozmer_heroy*-1}
            if(rozmer_heroz<0){rozmer_heroz=rozmer_heroz*-1} 
            //yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
            
            var rozmerx = rozmer_herox + rozmer_groundx;
            var rozmery = rozmer_heroy + rozmer_groundy;
            var rozmerz = rozmer_heroz + rozmer_groundz;
            
            
            if ((posx <= rozmerx && posy <= rozmery && posz <= rozmerz))        // pokud jsou vzdalenosti mezi objekty mensi nez soucet
            {                                                                   // jejich rozmeru v dane ose jsou objekty blizko
                return true;                                                    // dost na to aby se kolize testovala presnejším algoritmem
            } 
            else
            {
                return false;
            }
    }
    
    //studijní materiály: [1][3][5][6][7]

    function rotace_pro_5_alg(obj1, obj1_hrany, obj2, obj2_hrany)             //funkce ktera zrotuje objekty a otestuje je na kolizi
    {                                                                         // využívá ji algorimtus s OBB
        var euler_obj2 = obj2.rotationQuaternion.toEulerAngles();             //autor: Tomáš Pribyl
            
            for (let n=0; n < obj1_hrany.length; n++)                         // nejdríve se body posunou aby se rotace provádela podle správného body
            {
                obj1_hrany[n][0] = obj1_hrany[n][0] - obj1.position.x;
                obj1_hrany[n][1] = obj1_hrany[n][1] - obj1.position.y;
                obj1_hrany[n][2] = obj1_hrany[n][2] - obj1.position.z;
            } 
            
            var euler_obj1 = obj1.rotationQuaternion.toEulerAngles();            // poté se provede rotace ve všech trech osách
            
            var obj1_sin_x = Math.sin(euler_obj1.x);
            var obj1_cos_x = Math.cos(euler_obj1.x);
    
            for (let n=0; n < obj1_hrany.length; n++) 
            {
                let y = obj1_hrany[n][1];
                let z = obj1_hrany[n][2];
                obj1_hrany[n][1] = y * obj1_cos_x - z * obj1_sin_x;
                obj1_hrany[n][2] = z * obj1_cos_x + y * obj1_sin_x;
            }
            
            var obj1_sin_y = Math.sin(euler_obj1.y);
            var obj1_cos_y = Math.cos(euler_obj1.y);
    
            for (let n=0; n < obj1_hrany.length; n++) 
            {
                let x = obj1_hrany[n][0];
                var z = obj1_hrany[n][2];
                obj1_hrany[n][0] = x * obj1_cos_y + z * obj1_sin_y;
                obj1_hrany[n][2] = z * obj1_cos_y - x * obj1_sin_y;
            }
            
            var obj1_sin_z = Math.sin(euler_obj1.z);
            var obj1_cos_z = Math.cos(euler_obj1.z);  
    
            for (let n=0; n < obj1_hrany.length; n++) 
            {
                let x = obj1_hrany[n][0];
                let y = obj1_hrany[n][1];
                obj1_hrany[n][0] = x * obj1_cos_z - y * obj1_sin_z;
                obj1_hrany[n][1] = y * obj1_cos_z + x * obj1_sin_z;
            }  
            
            for (let n=0; n < obj1_hrany.length; n++)                              // body se vrati do puvodni pozice
            {
                obj1_hrany[n][0] = obj1_hrany[n][0] + obj1.position.x;
                obj1_hrany[n][1] = obj1_hrany[n][1] + obj1.position.y;
                obj1_hrany[n][2] = obj1_hrany[n][2] + obj1.position.z;
            } 
            
            
            for (let n=0; n < obj1_hrany.length; n++)                            // body se daji opet do nove pozice kvuli rotaci podle obj2 objektu
            {
                obj1_hrany[n][0] = obj1_hrany[n][0] - obj2.position.x;
                obj1_hrany[n][1] = obj1_hrany[n][1] - obj2.position.y;
                obj1_hrany[n][2] = obj1_hrany[n][2] - obj2.position.z;
            } 
            
            var obj1_sin_x = Math.sin(-euler_obj2.x);                          // opet se provedou rotace ve vsech trech osách
            var obj1_cos_x = Math.cos(-euler_obj2.x);                        
    
            for (let n=0; n < obj1_hrany.length; n++) 
            {
                let y = obj1_hrany[n][1];
                let z = obj1_hrany[n][2];
                obj1_hrany[n][1] = y * obj1_cos_x - z * obj1_sin_x;
                obj1_hrany[n][2] = z * obj1_cos_x + y * obj1_sin_x;
            }
            
            var obj1_sin_y = Math.sin(-euler_obj2.y);
            var obj1_cos_y = Math.cos(-euler_obj2.y);
    
            for (let n=0; n < obj1_hrany.length; n++) 
            {
                let x = obj1_hrany[n][0];
                var z = obj1_hrany[n][2];
                obj1_hrany[n][0] = x * obj1_cos_y + z * obj1_sin_y;
                obj1_hrany[n][2] = z * obj1_cos_y - x * obj1_sin_y;
            }
            
            var obj1_sin_z = Math.sin(-euler_obj2.z);
            var obj1_cos_z = Math.cos(-euler_obj2.z);  
    
            for (let n=0; n < obj1_hrany.length; n++) 
            {
                let x = obj1_hrany[n][0];
                let y = obj1_hrany[n][1];
                obj1_hrany[n][0] = x * obj1_cos_z - y * obj1_sin_z;
                obj1_hrany[n][1] = y * obj1_cos_z + x * obj1_sin_z;
            }
            
            for (let n=0; n < obj1_hrany.length; n++)                            // a body se opet vrátí do normálu
            {
                obj1_hrany[n][0] = obj1_hrany[n][0] + obj2.position.x;
                obj1_hrany[n][1] = obj1_hrany[n][1] + obj2.position.y;
                obj1_hrany[n][2] = obj1_hrany[n][2] + obj2.position.z;
            }
            

            for (let n=0; n < obj1_hrany.length; n++)      // pokud se nejaky bod obj1 nachazí uvnitr prostoru obj2 je to kolize
            {
                if(obj1_hrany[n][1] >= obj2_hrany[0][1] && obj1_hrany[n][1] <= obj2_hrany[2][1] && obj1_hrany[n][2] >= obj2_hrany[0][2] && obj1_hrany[n][2] <= obj2_hrany[1][2] )
                {
                    if(obj1_hrany[n][0] >= obj2_hrany[2][0] && obj1_hrany[n][0] <= obj2_hrany[6][0] && obj1_hrany[n][2] >= obj2_hrany[2][2] && obj1_hrany[n][2] <= obj2_hrany[3][2] )
                    {
                        if(obj1_hrany[n][0] >= obj2_hrany[0][0] && obj1_hrany[n][0] <= obj2_hrany[6][0] && obj1_hrany[n][1] >= obj2_hrany[0][1] && obj1_hrany[n][1] <= obj2_hrany[2][1] )
                        {
                            return 1;
                        }
                    }
                
                }
            }
            return 0;
        
    } 
    
    /*    ****************************************************************************************************
tato funkce doTrianglesIntersect je prevzatá od tretí strany
je dostupná na adrese github.com/kenny-evitt/three-js-triangle-triangle-collision-detection/blob/master/collision-tests.js
slouží detekci kolizí mezi trojúhleníky ve 3D prostoru pomocí separating axis theoremu
*/
    
    function doTrianglesIntersect(t1, t2) {

        var A0 = t1[0].clone();       // nejdríve se získají potrebné údaje -
        var A1 = t1[1].clone();       // body a vektory trojúhelníku a zapíší se do promenných
        var A2 = t1[2].clone();

        var E0 = A1.clone().subtract(A0);
        var E1 = A2.clone().subtract(A0);

        var E2 = E1.clone().subtract(E0);

        var N = E0.clone().cross(E1);

        var B0 = t2[0].clone();
        var B1 = t2[1].clone();
        var B2 = t2[2].clone();

        var F0 = B1.clone().subtract(B0);
        var F1 = B2.clone().subtract(B0);

        var F2 = F1.clone().subtract(F0);

        var M = F0.clone().cross(F1);


        var D = B0.clone().subtract(A0);


        function areProjectionsSeparated(p0, p1, p2, q0, q1, q2) {    // funkce zkoumá zda projekce na osách se protínají
            var min_p = Math.min(p0, p1, p2),
            max_p = Math.max(p0, p1, p2),
            min_q = Math.min(q0, q1, q2),
            max_q = Math.max(q0, q1, q2);

            return ((min_p > max_q) || (max_p < min_q));   // pokud se hodnoty mezi body neprotínají pak trojúhelníky nejsou v kolizi
        }

        {
            var p0 = 0,         // nyní se vytvorí celkem 11 os a vytvorí se na nich projekce bodu trojúhelníku
            p1 = 0,
            p2 = 0,
            q0 = N.dot(D),
            q1 = q0 + N.dot(F0),
            q2 = q0 + N.dot(F1);

            if (areProjectionsSeparated(p0, p1, p2, q0, q1, q2))     // tyto projekce se poté otestují zda se protínají
                return false;            // pokud se body neprotínají i jen na jedné ose, algoritmus koncí a vrací false, protože trojúhelníky v kolizi nejsou 
        }

        {
            var p0 = 0,
            p1 = M.dot(E0),
            p2 = M.dot(E1),
            q0 = M.dot(D),
            q1 = q0,
            q2 = q0;

            if (areProjectionsSeparated(p0, p1, p2, q0, q1, q2))
                return false;
        }

        {
            var p0 = 0,
            p1 = 0,
            p2 = -(N.dot(F0)),
            q0 = E0.clone().cross(F0).dot(D),
            q1 = q0,
            q2 = q0 + M.dot(E0);

            if (areProjectionsSeparated(p0, p1, p2, q0, q1, q2))
                return false;
        }

        {
            var p0 = 0,
            p1 = 0,
            p2 = -(N.dot(F1)),
            q0 = E0.clone().cross(F1).dot(D),
            q1 = q0 - M.dot(E0),
            q2 = q0;

            if (areProjectionsSeparated(p0, p1, p2, q0, q1, q2))
                return false;
        }

        {
            var p0 = 0,
            p1 = 0,
            p2 = -(N.dot(F2)),
            q0 = E0.clone().cross(F2).dot(D),
            q1 = q0 - M.dot(E0),
            q2 = q1;

            if (areProjectionsSeparated(p0, p1, p2, q0, q1, q2))
                return false;
        }

        {
            var p0 = 0,
            p1 = N.dot(F0),
            p2 = 0,
            q0 = E1.clone().cross(F0).dot(D),
            q1 = q0,
            q2 = q0 + M.dot(E1);

            if (areProjectionsSeparated(p0, p1, p2, q0, q1, q2))
                return false;
        }

        {
            var p0 = 0,
            p1 = N.dot(F1),
            p2 = 0,
            q0 = E1.clone().cross(F1).dot(D),
            q1 = q0 - M.dot(E1),
            q2 = q0;

            if (areProjectionsSeparated(p0, p1, p2, q0, q1, q2))
                return false;
        }

        {
            var p0 = 0,
            p1 = N.dot(F2),
            p2 = 0,
            q0 = E1.clone().cross(F2).dot(D),
            q1 = q0 - M.dot(E1),
            q2 = q1;

            if (areProjectionsSeparated(p0, p1, p2, q0, q1, q2))
                return false;
        }

        {
            var p0 = 0,
            p1 = N.dot(F0),
            p2 = p1,
            q0 = E2.clone().cross(F0).dot(D),
            q1 = q0,
            q2 = q0 + M.dot(E2);

            if (areProjectionsSeparated(p0, p1, p2, q0, q1, q2))
                return false;
        }

        {
            var p0 = 0,
            p1 = N.dot(F1),
            p2 = p1,
            q0 = E2.clone().cross(F1).dot(D),
            q1 = q0 - M.dot(E2),
            q2 = q0;

            if (areProjectionsSeparated(p0, p1, p2, q0, q1, q2))
                return false;
        }

        {
            var p0 = 0,
            p1 = N.dot(F2),
            p2 = p1,
            q0 = E2.clone().cross(F2).dot(D),
            q1 = q0 - M.dot(E2),
            q2 = q1;

            if (areProjectionsSeparated(p0, p1, p2, q0, q1, q2))
                return false;
        }

        return true;       // pokud se projekce na všech osách protínají vrací true a trojúhelníky tak v kolizi jsou
    }    // konec prevzate funkce
    
    //studijní materiály: [1][2][6][7]
    
    function doTrianglesIntersect2(t1, t2) {      // implementována funkce pro detekci kolizí mezi trojúhelníky
                                                  // autor: Tomáš Pribyl
    var A1 = t1[0].clone();
    var A2 = t1[1].clone();
    var A3 = t1[2].clone();
  
    var B1 = t2[0].clone();
    var B2 = t2[1].clone();
    var B3 = t2[2].clone();
 
  
    var trojuhelnik1 = [];                         // hodnoty se dají do polí pro lepší manipulaci                
    trojuhelnik1.push([A1.x, A1.y, A1.z]);                
    trojuhelnik1.push([A2.x, A2.y, A2.z]);                
    trojuhelnik1.push([A3.x, A3.y, A3.z]);                
    var trojuhelnik2 = [];                                         
    trojuhelnik2.push([B1.x, B1.y, B1.z]);             
    trojuhelnik2.push([B2.x, B2.y, B2.z]);           
    trojuhelnik2.push([B3.x, B3.y, B3.z]);  
    
    var trojuhelnik11 = [];                         // pole jejichz hodnoty se budou menit                
    trojuhelnik11.push([A1.x, A1.y, A1.z]);          // a do kterych se zapise rotace      
    trojuhelnik11.push([A2.x, A2.y, A2.z]);                
    trojuhelnik11.push([A3.x, A3.y, A3.z]);                
    var trojuhelnik22 = [];                                         
    trojuhelnik22.push([B1.x, B1.y, B1.z]);             
    trojuhelnik22.push([B2.x, B2.y, B2.z]);           
    trojuhelnik22.push([B3.x, B3.y, B3.z]);                 

     var xx,yy;
      
    
    
for (var index1 = 10; index1<=180; index1 = index1+10){         // dvojity for
for (var index2 = 10; index2<=180; index2 = index2+10){         // kazdy pro rotaci v jine ose
  
    
    var uhelsin1 = Math.sin(index1*(Math.PI/180));           // uhly
    var uhelcos1 = Math.cos(index1*(Math.PI/180)); 
    
    var uhelsin2 = Math.sin(index2*(Math.PI/180));
    var uhelcos2 = Math.cos(index2*(Math.PI/180)); 
    
    for(var x = 0; x<3; x++)
    {
        xx = trojuhelnik1[x][0];                                    // rotace 1. trojuhelniku v ose z
        yy = trojuhelnik1[x][1];
        trojuhelnik11[x][0] = xx * uhelcos1 - yy * uhelsin1;    
        trojuhelnik11[x][1] = yy * uhelcos1 + xx * uhelsin1;
        trojuhelnik11[x][2] = trojuhelnik1[x][2]; 
        
        xx = trojuhelnik11[x][1];                                  // rotace 1. trojuhelniku v ose x
        yy = trojuhelnik11[x][2];
        trojuhelnik11[x][0] = trojuhelnik11[x][0] ;    
        trojuhelnik11[x][1] = xx * uhelcos2 - yy * uhelsin2; 
        trojuhelnik11[x][2] = yy * uhelcos2 + xx * uhelsin2; 

        xx = trojuhelnik2[x][0];                                    // rotace 2. trojuhelniku v ose z
        yy = trojuhelnik2[x][1];
        trojuhelnik22[x][0] = xx * uhelcos1 - yy * uhelsin1;    
        trojuhelnik22[x][1] = yy * uhelcos1 + xx * uhelsin1;
        trojuhelnik22[x][2] = trojuhelnik2[x][2]; 
       
        xx = trojuhelnik22[x][1];                                       // rotace 2. trojuhelniku v ose x
        yy = trojuhelnik22[x][2]; 
        trojuhelnik22[x][0] = trojuhelnik22[x][0] ;    
        trojuhelnik22[x][1] = xx * uhelcos2 - yy * uhelsin2; 
        trojuhelnik22[x][2] = yy * uhelcos2 + xx * uhelsin2;
        
    }
    
    var min_p = Math.min(trojuhelnik22[0][0], trojuhelnik22[1][0], trojuhelnik22[2][0]),     // kontrola zda jsou body 1. troj v ose x vyse 
    max_p = Math.max(trojuhelnik22[0][0], trojuhelnik22[1][0], trojuhelnik22[2][0]),        // nebo nize nez body 2. troj
    min_q = Math.min(trojuhelnik11[0][0], trojuhelnik11[1][0], trojuhelnik11[2][0]),
    max_q = Math.max(trojuhelnik11[0][0], trojuhelnik11[1][0], trojuhelnik11[2][0]);
    
    if ((min_p > max_q) || (max_p < min_q))
    {
        return false                                                                       // pokud ano, tak trojúhelníky nejsou v kolizi
    }
    
    min_p = Math.min(trojuhelnik22[0][1], trojuhelnik22[1][1], trojuhelnik22[2][1]),      // to, same ale osa y
    max_p = Math.max(trojuhelnik22[0][1], trojuhelnik22[1][1], trojuhelnik22[2][1]),
    min_q = Math.min(trojuhelnik11[0][1], trojuhelnik11[1][1], trojuhelnik11[2][1]),
    max_q = Math.max(trojuhelnik11[0][1], trojuhelnik11[1][1], trojuhelnik11[2][1]);
    
    if ((min_p > max_q) || (max_p < min_q))
    {
        return false
    }
    
    min_p = Math.min(trojuhelnik22[0][2], trojuhelnik22[1][2], trojuhelnik22[2][2]),     // to, same ale osa z
    max_p = Math.max(trojuhelnik22[0][2], trojuhelnik22[1][2], trojuhelnik22[2][2]),
    min_q = Math.min(trojuhelnik11[0][2], trojuhelnik11[1][2], trojuhelnik11[2][2]),
    max_q = Math.max(trojuhelnik11[0][2], trojuhelnik11[1][2], trojuhelnik11[2][2]);
    
    if ((min_p > max_q) || (max_p < min_q))
    {
        return false
    }
    
    }}
    
    return true;                              // jinak jsou trojuhelníky v kolizi
  
}
    
    BABYLON.Vector3.prototype.dot = function (b)        // funkce pro výpocet vektoru
     {
        return this.x * b.x + this.y * b.y + this.z * b.z;
     }