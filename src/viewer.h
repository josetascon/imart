#ifndef __VIEWER_HPP__
#define __VIEWER_HPP__

// VTK Libraries
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkImageFlip.h>
#include <vtkImageImport.h>
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
// #include <vtkCamera.h>
// #include <vtkProperty.h>

// Local headers
#include "object.h"

namespace imart
{

template <typename type>
class viewer : public inherit<viewer<type>, object>
{
public:
    //Type definitions
    using self    = viewer;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    int num_plots;
    int added_plots;
    std::vector<bool>                   _flips_;
    std::vector<bool>                   _update_;
    std::vector<unsigned int>           _size_;
    std::vector<unsigned int>           _offset_;
    std::vector<unsigned int>           _subplot_;
    std::vector<std::vector<double>>    _viewports_;
    std::vector<std::string>            _masks_;
    std::vector<std::string>            _titles_;

    std::vector<std::shared_ptr<type>>              _images_;
    std::vector<vtkSmartPointer<vtkImageData>>      _vtkimg_;
    std::vector<vtkSmartPointer<vtkImageActor>>     _actors_;
    std::vector<vtkSmartPointer<vtkRenderer>>       _renders_;
    vtkSmartPointer<vtkRenderWindow>                _renwin_;
    vtkSmartPointer<vtkRenderWindowInteractor>      _interactor_;
    vtkSmartPointer<vtkInteractorStyleImage>        _style_;

    // // std::vector<vtkVolume *>    _volumes_;
    // std::vector< vtkSmartPointer<vtkRenderer> >                 _renders_;
    // std::vector< vtkSmartPointer<vtkRenderWindow> >             _renwins_;
    // std::vector< vtkSmartPointer<vtkRenderWindowInteractor> >   _interactors_;

    // ===========================================
    // Functions
    // ===========================================
    vtkSmartPointer<vtkImageData> imart2vtk(std::shared_ptr<type> image_pointer);

public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    // Constructors
    //Constructor
    viewer();
    //Destructor
    ~viewer();

    // ===========================================
    // Set Functions
    // ===========================================
    // Change window size
    void size(unsigned int h, unsigned int v);

    // ===========================================
    // Functions
    // ===========================================
    // Define subplots in window
    void subplot(unsigned int row, unsigned int col);
    // Add image to window
    void add_image(std::shared_ptr<type> image_pointer);
    void add_image(std::shared_ptr<type> image_pointer, std::string title, std::string mask, bool flip);
    // Visualize all the image added
    void setup();
    void update_image(std::shared_ptr<type> img, int id);
    // void update(images);

    void visualize();
    void show();
};


// ===========================================
//      Functions of Class registration
// ===========================================

// ===========================================
// Constructor Functions
// ===========================================
// Constructor
template <typename type>
viewer<type>::viewer()
{
    num_plots = 1;
    added_plots = 0;
    // _size_ = std::vector<unsigned int>{1200, 900};
    _size_ = std::vector<unsigned int>{700, 450};
    _offset_ = std::vector<unsigned int>{100, 100};
    _subplot_ = std::vector<unsigned int>{1, 1};
    _renwin_ = vtkSmartPointer<vtkRenderWindow>::New();
    _style_ = vtkSmartPointer<vtkInteractorStyleImage>::New();
    _interactor_ = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    _renwin_->SetSize(_size_[0], _size_[1]);
};

// Declaration of method: viewer::~viewer()
template <typename type>
viewer<type>::~viewer() 
{
    ;
};

// Declaration of method: viewer::size()
template <typename type>
void viewer<type>::size(unsigned int h, unsigned int v)
{
    _size_[0] = h;
    _size_[1] = v;
    _renwin_->SetSize(_size_[0], _size_[1]); 
};

// ===========================================
// Functions
// ===========================================
// Declaration of method subplot
template <typename type>
void viewer<type>::subplot(unsigned int row, unsigned int col)
{
    _subplot_[0] = row;
    _subplot_[1] = col;
    num_plots = row*col;
    //std::cout << "Row: " << _subplot_[0] << ", Col: " << _subplot_[1] << std::endl;
};

// Declaration of method: void viewer::add_image(type *p_image)
template <typename type>
void viewer<type>::add_image(std::shared_ptr<type> image_pointer)
{
    std::stringstream title;
    title << "Image" << (added_plots + 1);
    std::string mask = "Grey";
    add_image(image_pointer, title.str(), mask, true);
};

// Declaration of method: viewer::add_image(type *p_image, bool flip, std::string title, std::string mask)
template <typename type>
void viewer<type>::add_image(std::shared_ptr<type> image_pointer, std::string title, std::string mask, bool flip)
{
    _images_.push_back(image_pointer);
    _flips_.push_back(flip);
    _titles_.push_back(title);
    _masks_.push_back(mask);
    added_plots++;
};

template <typename type>
vtkSmartPointer<vtkImageData> viewer<type>::imart2vtk(std::shared_ptr<type> image_pointer)
{
    vtkSmartPointer<vtkImageImport> imageImport = vtkSmartPointer<vtkImageImport>::New();
    auto dim = image_pointer->get_dimension();
    auto size = image_pointer->get_size();
    auto spacing = image_pointer->get_spacing();
    auto origin = image_pointer->get_origin();

    if (dim == 2)
    {
        imageImport->SetDataSpacing(spacing[0], spacing[1], 1);
        imageImport->SetDataOrigin(origin[0], origin[1], 0);
        imageImport->SetWholeExtent(0, size[0]-1, 0, size[1]-1, 0, 0);
    }
    else if (dim == 3)
    {
        imageImport->SetDataSpacing(spacing[0], spacing[1], spacing[2]);
        imageImport->SetDataOrigin(origin[0], origin[1], origin[2]);
        imageImport->SetWholeExtent(0, size[0]-1, 0, size[1]-1, 0, size[2]-1);
    };
    imageImport->SetDataExtentToWholeExtent();
    // TODO: modify this according to type
    imageImport->SetDataScalarTypeToDouble();
    // if (strcmp(image_pointer->get_type().c_str(),"unsigned_char"))
    //     imageImport->SetDataScalarTypeToUnsignedChar();
    // else if (strcmp(image_pointer->get_type().c_str(),"unsigned_short"))
    //     imageImport->SetDataScalarTypeToUnsignedShort();
    // else if (strcmp(image_pointer->get_type().c_str(),"short"))
    //     imageImport->SetDataScalarTypeToShort();
    // else if (strcmp(image_pointer->get_type().c_str(),"float"))
    //     imageImport->SetDataScalarTypeToFloat();
    // else if (strcmp(image_pointer->get_type().c_str(),"double"))
    //     imageImport->SetDataScalarTypeToDouble();
    // imageImport->SetNumberOfScalarComponents(image_pointer->get_channels());
    
    // std::vector<double> vec = image_pointer->get_data()->std_vector();
    // imageImport->SetImportVoidPointer(vec.data());

    imageImport->SetImportVoidPointer(image_pointer->get_data()->data());    
    imageImport->Update();

    // return imageImport->GetOutput();

    vtkSmartPointer<vtkImageFlip> flip = vtkSmartPointer<vtkImageFlip>::New();
    flip->SetInputData(imageImport->GetOutput());
    flip->SetFilteredAxis(1);
    flip->Update();
    return flip->GetOutput();
};

// Declaration of method: void viewer::visualize()
template <typename type>
void viewer<type>::setup()
{
    if (num_plots < added_plots)
    {
        throw "Error. Added more plots than the specified";
    }

    unsigned int len = _images_.size();
    unsigned int row = _subplot_[0];
    unsigned int col = _subplot_[1];
    unsigned int win_h = _size_[0];
    unsigned int win_v = _size_[1];

    std::vector<std::vector<unsigned int>> positions;
    for (unsigned int i = 0; i < row; i++)
    {
        for (unsigned int j = 0; j < col; j++)
        {
            unsigned int x = _offset_[1]+j*win_h;
            unsigned int y = _offset_[0]+i*win_v;
            positions.push_back(std::vector<unsigned int>{x,y});
        };
    };

    // Compute the viewports based on subplot
    for (double i = 0.0; i < row; i++)
    {
        for (double j = 0.0; j < col; j++)
        {
            // Viewport format {xmin,ymin,xmax,ymax}
            // std::cout << "viewport " << std::endl;
            // std::cout << static_cast<double>(j/col) << " " << static_cast<double>(i/row) << " " << 
            //      static_cast<double>((j+1.0)/col) << " " << static_cast<double>((i+1.0)/row) << "\n";
            std::vector<double> viewport = {static_cast<double>(j/col), static_cast<double> (i/row),
                                 static_cast<double> ((j+1.0)/col), static_cast<double> ((i+1.0)/row)};
            _viewports_.push_back(viewport);
        };
    };

    // printf("loop imart2vtk\n");
    for (unsigned int k = 0; k < len ; k++)
    {
        // Get vtk image from image data
        // printf("image %d\n", k);
        vtkSmartPointer<vtkImageData> img = imart2vtk(_images_[k]);
        _vtkimg_.push_back(img);
        // _images_[k]->print();

        //Actor
        vtkSmartPointer<vtkImageActor> actor = vtkSmartPointer<vtkImageActor>::New();
        _actors_.push_back(actor);
        actor->GetMapper()->SetInputData(img);

        // Renderer
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
        _renders_.push_back(renderer);
        renderer->SetViewport(_viewports_[k].data());
        renderer->AddActor(actor);
        // renderer->SetBackground(0.3,0.5,0.8); // blue
        renderer->SetBackground(0,0,0); // dark
        renderer->ResetCamera();
        _renwin_->AddRenderer(renderer);
    };

    return;
};

template <typename type>
void viewer<type>::update_image(std::shared_ptr<type> img, int id)
{
    if (id < _images_.size())
    {
        _images_[id] = img;
        _vtkimg_[id] = imart2vtk(_images_[id]);
        _actors_[id]->GetMapper()->SetInputData(_vtkimg_[id]);
        // _renders_[id]->SetViewport(_viewports_[id].data());
        // _renders_[id]->AddActor(_actors_[id]);
        // _renders_[id]->SetBackground(0,0,0); // dark
        // _renders_[id]->ResetCamera();
        // _renwin_->AddRenderer(_renders_[id]);
    }
    else
    {
        printf("[Warning][Viewer] Invalid update due to image id");
    };
    // visualize();
};

template <typename type>
void viewer<type>::visualize()
{
    _renwin_->Render();
    // _interactor_->SetRenderWindow(_renwin_);
    // _interactor_->Initialize();
    // _interactor_->Start();
};

template <typename type>
void viewer<type>::show()
{
    printf("Running Viewer\n");
    _interactor_->SetRenderWindow(_renwin_);
    _interactor_->Initialize();
    _interactor_->Start();
};



}; //end namespace

#endif