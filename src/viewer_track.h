#ifndef __viewer_track_TRACK_HPP__
#define __viewer_track_TRACK_HPP__

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
// #include <vtkProperty.h>

#include <vtkNamedColors.h>
#include <vtkImageMapToColors.h>
#include <vtkLookupTable.h>

#include <vtkCamera.h>


// Local headers
#include "object.h"

namespace imart
{

template <typename type>
class viewer_track : public inherit<viewer_track<type>, object>
{
public:
    //Type definitions
    using self    = viewer_track;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    int added_images;
    std::vector<unsigned int>           _size_;
    std::vector<unsigned int>           _offset_;
    std::vector<std::string>            _color_names_;

    std::vector<std::shared_ptr<type>>              _images_;
    std::vector<vtkSmartPointer<vtkImageData>>      _vtkimages_;
    std::vector<vtkSmartPointer<vtkImageActor>>     _actors_;
    vtkSmartPointer<vtkNamedColors>                 _colors_;
    vtkSmartPointer<vtkRenderer>                    _render_;
    vtkSmartPointer<vtkRenderWindow>                _renwin_;
    vtkSmartPointer<vtkRenderWindowInteractor>      _interactor_;
    vtkSmartPointer<vtkInteractorStyleImage>        _style_;

    // ===========================================
    // Functions
    // ===========================================
    void init_colors();
    vtkSmartPointer<vtkImageData> imart2vtk(std::shared_ptr<type> image_pointer);

public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    // Constructors
    //Constructor
    viewer_track();
    //Destructor
    ~viewer_track();

    // ===========================================
    // Set Functions
    // ===========================================
    // Change window size
    void size(unsigned int h, unsigned int v);

    // ===========================================
    // Functions
    // ===========================================
    // Add image
    void add_image(std::shared_ptr<type> image_pointer);
    // Update image data
    void update(int id);
    // Visualize all the images added
    void setup();
    void render();
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
viewer_track<type>::viewer_track()
{
    added_images = 0;
    // _size_ = std::vector<unsigned int>{1200, 900};
    _size_ = std::vector<unsigned int>{700, 450};
    _offset_ = std::vector<unsigned int>{100, 100};
    _colors_ = vtkSmartPointer<vtkNamedColors>::New();
    _render_ = vtkSmartPointer<vtkRenderer>::New();
    _renwin_ = vtkSmartPointer<vtkRenderWindow>::New();
    _style_ = vtkSmartPointer<vtkInteractorStyleImage>::New();
    _interactor_ = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    _renwin_->SetSize(_size_[0], _size_[1]);

    init_colors();
};

// Declaration of method: viewer_track::~viewer_track()
template <typename type>
viewer_track<type>::~viewer_track() 
{
    ;
};

template <typename type>
void viewer_track<type>::init_colors()
{
    _color_names_ = {"Crimson", "LimeGreen", "RoyalBlue", "Khaki", "Violet", "Sienna"};
};

// Declaration of method: viewer_track::size()
template <typename type>
void viewer_track<type>::size(unsigned int h, unsigned int v)
{
    _size_[0] = h;
    _size_[1] = v;
    _renwin_->SetSize(_size_[0], _size_[1]); 
};

// ===========================================
// Functions
// ===========================================

// Declaration of method: void viewer_track::add_image(type *p_image)
template <typename type>
void viewer_track<type>::add_image(std::shared_ptr<type> image_pointer)
{
    _images_.push_back(image_pointer);
    added_images++;
};


template <typename type>
// vtkSmartPointer<vtkImageImport> viewer_track<type>::imart2vtk(std::shared_ptr<type> image_pointer)
vtkSmartPointer<vtkImageData> viewer_track<type>::imart2vtk(std::shared_ptr<type> image_pointer)
{   
    image_pointer->to_cpu(); // update in case gpu image
    vtkSmartPointer<vtkImageImport> imageImport = vtkSmartPointer<vtkImageImport>::New();
    auto dim = image_pointer->get_dimension();
    auto size = image_pointer->get_size();
    auto spacing = image_pointer->get_spacing();
    auto origin = image_pointer->get_origin();

    // double direction[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    // imageImport->SetDataDirection(direction);

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
    imageImport->SetDataScalarTypeToFloat();
    // imageImport->SetDataScalarTypeToDouble();
    // imageImport->SetDataScalarTypeToUnsignedChar();


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

    return imageImport->GetOutput();

    // vtkSmartPointer<vtkImageFlip> flip = vtkSmartPointer<vtkImageFlip>::New();
    // flip->SetInputData(imageImport->GetOutput());
    // flip->SetFilteredAxis(1);
    // flip->Update();
    // return flip->GetOutput();
    // // return imageImport->GetOutput();
    // return imageImport;
};

// Declaration of method: void viewer_track::visualize()
template <typename type>
void viewer_track<type>::setup()
{
    unsigned int len = _images_.size();
    
    // printf("loop imart2vtk\n");
    for (unsigned int k = 0; k < len ; k++)
    {
        // Get vtk image from image data
        vtkSmartPointer<vtkImageData> img = imart2vtk(_images_[k]);
        _vtkimages_.push_back(img);

        //Actor
        vtkSmartPointer<vtkImageActor> actor = vtkSmartPointer<vtkImageActor>::New();
        _actors_.push_back(actor);

        if (k == 0)
        {
            vtkSmartPointer<vtkImageFlip> flip = vtkSmartPointer<vtkImageFlip>::New();
            flip->SetInputData(img);
            flip->SetFilteredAxis(1);
            flip->Update();
            actor->GetMapper()->SetInputConnection(flip->GetOutputPort());
            // actor->SetInputData(img);
        }
        else
        {
            auto dc = _colors_->GetColor3d(_color_names_[k-1]);
            // std::cout << _color_names_[k-1] << std::endl;
            // printf("Color: %4.3f, %4.3f, %4.3f\n", dc[0], dc[1], dc[2]);
            
            vtkSmartPointer<vtkLookupTable> lookup_table =
                vtkSmartPointer<vtkLookupTable>::New();
            lookup_table->SetNumberOfTableValues(2);
            lookup_table->SetRange(0.0, 1.0);
            lookup_table->SetTableValue( 0, 0.0, 0.0, 0.0, 0.0 ); //label 0 is transparent
            lookup_table->SetTableValue( 1, dc[0], dc[1], dc[2], 0.3 ); //label 1 with color
            lookup_table->Build();

            // Pass the original image and the lookup table to a filter to create a color image:
            vtkSmartPointer<vtkImageMapToColors> mask_to_color =
                vtkSmartPointer<vtkImageMapToColors>::New();
            mask_to_color->SetLookupTable(lookup_table);
            mask_to_color->PassAlphaToOutputOn();
            mask_to_color->SetInputData( img );
            mask_to_color->Update();

            vtkSmartPointer<vtkImageFlip> flip = vtkSmartPointer<vtkImageFlip>::New();
            flip->SetInputConnection( mask_to_color->GetOutputPort() );
            flip->SetFilteredAxis(1);
            flip->Update();
            actor->GetMapper()->SetInputConnection(flip->GetOutputPort());

            actor->GetMapper()->SetInputConnection( flip->GetOutputPort() );
            // actor->GetMapper()->SetInputConnection( mask_to_color->GetOutputPort() );
        }

        // Renderer
        _render_->AddActor(actor);
    };

    // Renderer
    _render_->SetBackground(0,0,0); // dark
    // _render_->SetBackground(0.3,0.5,0.8); // blue
    _render_->ResetCamera();
    _renwin_->AddRenderer(_render_);

    

    return;
};

template <typename type>
void viewer_track<type>::update( int id )
{
    if (id < _images_.size())
        _vtkimages_[id]->Modified();
    else
        printf("[Warning][viewer_track] Invalid update due to invalid image id");
};

template <typename type>
void viewer_track<type>::render()
{
    _renwin_->Render();
};

template <typename type>
void viewer_track<type>::show()
{
    printf("Running viewer_track\n");
    _interactor_->SetRenderWindow(_renwin_);
    _interactor_->Initialize();
    _interactor_->Start();
};

}; //end namespace

#endif