#include "mainwindow.h"
#include "ui_mainwindow.h"


/* **** début de la partie boutons et IHM **** */


// exemple pour charger un fichier .obj
void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}

// exemple pour construire un mesh face par face
void MainWindow::on_pushButton_generer_clicked()
{
    test();
    /*MyMesh mesh;

    // on construit une liste de sommets
    MyMesh::VertexHandle sommets[8];
    sommets[0] = mesh.add_vertex(MyMesh::Point(-1, -1,  1));
    sommets[1] = mesh.add_vertex(MyMesh::Point( 1, -1,  1));
    sommets[2] = mesh.add_vertex(MyMesh::Point( 1,  1,  1));
    sommets[3] = mesh.add_vertex(MyMesh::Point(-1,  1,  1));
    sommets[4] = mesh.add_vertex(MyMesh::Point(-1, -1, -1));
    sommets[5] = mesh.add_vertex(MyMesh::Point( 1, -1, -1));
    sommets[6] = mesh.add_vertex(MyMesh::Point( 1,  1, -1));
    sommets[7] = mesh.add_vertex(MyMesh::Point(-1,  1, -1));


    // on construit des faces à partir des sommets

    std::vector<MyMesh::VertexHandle> uneNouvelleFace;

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[0]);
    uneNouvelleFace.push_back(sommets[1]);
    uneNouvelleFace.push_back(sommets[2]);
    uneNouvelleFace.push_back(sommets[3]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[7]);
    uneNouvelleFace.push_back(sommets[6]);
    uneNouvelleFace.push_back(sommets[5]);
    uneNouvelleFace.push_back(sommets[4]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[1]);
    uneNouvelleFace.push_back(sommets[0]);
    uneNouvelleFace.push_back(sommets[4]);
    uneNouvelleFace.push_back(sommets[5]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[2]);
    uneNouvelleFace.push_back(sommets[1]);
    uneNouvelleFace.push_back(sommets[5]);
    uneNouvelleFace.push_back(sommets[6]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[3]);
    uneNouvelleFace.push_back(sommets[2]);
    uneNouvelleFace.push_back(sommets[6]);
    uneNouvelleFace.push_back(sommets[7]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[0]);
    uneNouvelleFace.push_back(sommets[3]);
    uneNouvelleFace.push_back(sommets[7]);
    uneNouvelleFace.push_back(sommets[4]);
    mesh.add_face(uneNouvelleFace);

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);*/

}

/* **** fin de la partie boutons et IHM **** */



/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

void MainWindow::test()
{
    resetAllColorsAndThickness(&mesh);

    generateMatrix(&mesh, 10);
    /*  3x3
      y [   2   5   8   ]
      ^ [   1   4   7   ]
      | [   0   3   6   ]
        -> x
    */ //parametrisation et lissage

    printMatrix(&mesh);
    printVertices(&mesh);

    qDebug() << "nb vertices : : " << mesh.n_vertices() << "|| myMatrix size : : " << myMatrix.size();
    unsigned int id1 = 1;
    unsigned int id2 = 31;

    printVertex(&mesh, myMatrix.at(id1), MyMesh::Color(0, 255, 0)); // vert
    printVertex(&mesh, myMatrix.at(id2), MyMesh::Color(255, 0, 0)); // rouge
     qDebug() << "id matrix : : " << myMatrix.at(id1).idx();
    qDebug() << "p : : " << mesh.point(myMatrix.at(id1))[0] << mesh.point(myMatrix.takeAt(id1))[1];
    qDebug() << "p : : " << mesh.point(myMatrix.at(id2))[0] << mesh.point(myMatrix.takeAt(id2))[1];
    //qDebug() << "p : : " << mesh.point(myMatrix.takeAt(31))[0] << mesh.point(myMatrix.takeAt(31))[1];


    printVertex(&mesh, mesh.vertex_handle(id1), MyMesh::Color(255, 255, 0)); // jaune
    printVertex(&mesh, mesh.vertex_handle(id2), MyMesh::Color(0, 255, 255)); // cyan
    qDebug() << "p : : " << mesh.point(mesh.vertex_handle(id1))[0] << mesh.point(mesh.vertex_handle(id1))[1];
    qDebug() << "p : : " << mesh.point(mesh.vertex_handle(id2))[0] << mesh.point(mesh.vertex_handle(id2))[1];

    /*MyMesh::VertexHandle vh1 = mesh.vertex_handle(id1);
    MyMesh::VertexHandle vh2 = mesh.vertex_handle(id2);*/
    MyMesh::VertexHandle vh1 = myMatrix.at(id1);
    MyMesh::VertexHandle vh2 = myMatrix.at(id2);
    //lineBresenhamAlgorithm(&mesh, myMatrix.takeAt(0), myMatrix.takeAt(17));
    lineBresenhamAlgorithm(&mesh, vh1, vh2);
    displayMesh(&mesh);
}

void MainWindow::generateMatrix(MyMesh *_mesh, unsigned int square_matrix)
{
    /*matrix = Matrix(square_matrix);
    int id = 0;
    for (unsigned int x = 0; x < square_matrix; x++) {
        for (unsigned int y = 0; y < square_matrix; y++) {
            matrix.setTo(x, y, id);
            id ++;
        }
    }*/

    int i = 0;
    qDebug() << "myMatrix capacity : " << myMatrix.capacity();
    myMatrix.reserve( static_cast<int>( pow(square_matrix, 2) ) );
    for (float x = 0; x < square_matrix; x++) {
        for (float y = 0; y < square_matrix; y++) {
            //myMatrix.push_back(QVector3D(x, y, 0.0f));
            //myMatrix.push_back(_mesh->add_vertex(MyMesh::Point(x, y, 0.0f)));
            myMatrix.insert(i, _mesh->add_vertex(MyMesh::Point(x, y, 0.0f)));
            i++;
        }
    }
    qDebug() << "myMatrix capacity : " << myMatrix.capacity();
}

void MainWindow::printMatrix(MyMesh* _mesh)
{
    qDebug() << __FUNCTION__ << " :";
    QVectorIterator<MyMesh::VertexHandle> v_it(myMatrix);
    while (v_it.hasNext()) {
        MyMesh::VertexHandle vh = v_it.next();
        printVertex(_mesh, vh, MyMesh::Color(0, 0, 0));
        qDebug() << "id : : " << vh.idx() << " || x, y : : " << _mesh->point(vh)[0] << _mesh->point(vh)[1];
    }
}

void MainWindow::printVertices(MyMesh *_mesh)
{
    qDebug() << __FUNCTION__ << " :";
    for (MyMesh::VertexIter v_it = _mesh->vertices_sbegin(); v_it != _mesh->vertices_end(); v_it++) {
        /*_mesh->data(*v_it).thickness = 2;
        _mesh->set_color(*v_it, MyMesh::Color(0, 0, 0));*/
        printVertex(_mesh, *v_it, MyMesh::Color(0, 0, 0));
        qDebug() << "id : : " << v_it->idx() << " || x, y : : " << _mesh->point(*v_it)[0] << _mesh->point(*v_it)[1];
    }
}

void MainWindow::printVertex(MyMesh *_mesh, MyMesh::VertexHandle vh, MyMesh::Color color)
{
    _mesh->data(vh).thickness = 6;
    _mesh->set_color(vh, color);
}

void MainWindow::lineBresenhamAlgorithm(MyMesh* _mesh, MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2)
{
    /*  3x3
      y [   2   5   8   ]
      ^ [   1   4   7   ]
      | [   0   3   6   ]
        -> x
    */
    int x1 = static_cast<int>( _mesh->point(vh1)[0] );
    int y1 = static_cast<int>( _mesh->point(vh1)[1] );
    int x2 = static_cast<int>( _mesh->point(vh2)[0] );
    int y2 = static_cast<int>( _mesh->point(vh2)[1] );
    qDebug() << "x1, y1 : : " << x1 << y1 << " || x2, y2 : : " << x2 << y2;

    int x, y, dx, dy;
    double e, e_1_0, e_0_1; // error value and incrementals

    dy = static_cast<int>( y2 - y1 );
    dx = static_cast<int>( x2 - x1 );
    y = static_cast<int>( y1 ); // range init
    e = 0.0; // error value init
    e_1_0 = static_cast<double>( dy ) / static_cast<double> ( dx );
    e_0_1 = -1.0;
    qDebug() << "dx, dy : : " << dx << dy << dy / dx;
    qDebug() << "e : : " << e << " || e_1_0 : : " << e_1_0 << " || e_0_1 : : " << e_0_1;

    double size_matrix_square = sqrt(_mesh->n_vertices()/*myMatrix.capacity()*/);
    for (int x = x1; x <= x2; x++) {

        qDebug() << "x : : " << x << " || y : : " << y << " || e : : " << e << " || e_1_0 0_1 : : " << e_1_0 << e_0_1;
        unsigned int id_v = static_cast<unsigned int>( y + x * size_matrix_square );
        qDebug() << "id_v : : " << id_v;

        printVertex(_mesh, _mesh->vertex_handle(id_v), MyMesh::Color(0, 0, 255));

        // error for next vertex with even range
        e = e + e_1_0;
        qDebug() << "e : : " << e << " || e_1_0 : : " << e_1_0 << " || e_0_1 : : " << e_0_1;
        if ( e >= 0.5) {
            qDebug() << "------------------------";
            y = y + 1; // instead to choose the next vertex with even range
            e = e + e_0_1; // adjusting error in new range
        }
    }
    qDebug() << "--- end";
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }


    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}


