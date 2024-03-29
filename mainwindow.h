#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <QVector>
#include <QVector3D>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    EdgeAttributes( OpenMesh::Attributes::Color | OpenMesh::Attributes::Status );
    // vertex thickness
    VertexTraits{float thickness; float value;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1);
    void resetAllColorsAndThickness(MyMesh* _mesh);

    // --- test ---
    void test();

    // --- branch Function_generateMatrix ---
    void generateMatrix(MyMesh* _mesh, unsigned int square_matrix);
    void printMatrix(MyMesh* _mesh);
    void printVertices(MyMesh* _mesh);
    void printVertex(MyMesh* _mesh, MyMesh::VertexHandle vh, MyMesh::Color color);

private slots:
    void on_pushButton_chargement_clicked();

    void on_pushButton_generer_clicked();

private:

    MyMesh mesh;

    Ui::MainWindow *ui;

    //QVector<QVector3D> myMatrix;
    QVector<MyMesh::VertexHandle> myMatrix;
};

#endif // MAINWINDOW_H
