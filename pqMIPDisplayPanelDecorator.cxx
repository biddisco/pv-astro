/*=========================================================================

 Program:   Visualization Toolkit
 Module:    pqMIPDisplayPanelDecorator.cxx

 Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 All rights reserved.
 See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

 =========================================================================*/

// .NAME pqMIPDisplayPanelDecorator
// .SECTION Thanks
// <verbatim>
//
//  This file is part of the MIPs plugin developed and contributed by
//
// </verbatim>

#include "pqMIPDisplayPanelDecorator.h"
#include "ui_pqMIPDisplayPanelDecorator.h"

// Server Manager Includes.
#include "vtkCommand.h"
#include "vtkDataObject.h"
#include "vtkEventQtSlotConnect.h"
#include "vtkSmartPointer.h"
#include "vtkSMProperty.h"
#include "vtkSMPropertyHelper.h"
#include "vtkSMPVRepresentationProxy.h"
#include "vtkPVDataInformation.h"
#include "vtkSMStringVectorProperty.h"
#include "vtkSMEnumerationDomain.h"
#include "vtkPVDataSetAttributesInformation.h"
#include "vtkPVArrayInformation.h"

// Qt Includes.
#include <QVBoxLayout>
#include <QComboBox>
#include <QGroupBox>

// ParaView Includes.
#include "pqDisplayProxyEditor.h"
#include "pqRepresentation.h"
#include "pqFieldSelectionAdaptor.h"
#include "pqPropertyLinks.h"
#include "pqSMAdaptor.h"
#include "pqPipelineRepresentation.h"
#include "pqVariableType.h"
#include "pqScalarsToColors.h"
#include "pqWidgetRangeDomain.h"
#include "pqFieldSelectionAdaptor.h"
#include "pqSignalAdaptors.h"
#include "pqLookupTableManager.h"
#include "pqApplicationCore.h"
#include "pqCoreUtilities.h"

// enum ElementTypes{ INT, DOUBLE, STRING };

class pqMIPDisplayPanelDecorator::pqInternals: public Ui::pqMIPDisplayPanelDecorator
{
public:
  vtkSmartPointer<vtkEventQtSlotConnect> VTKConnect;
  pqPropertyLinks                        Links;
  vtkSMPVRepresentationProxy            *RepresentationProxy;
  pqPipelineRepresentation              *PipelineRepresentation;
  QWidget                               *Frame;
  QList<pqScalarsToColors*>              ColorTableList;
  int                                    TableIndex;
  //
  pqInternals(QWidget* parent)
  {
    this->VTKConnect     = vtkSmartPointer<vtkEventQtSlotConnect>::New();
    this->Frame          = 0;
    this->TableIndex     = 0;
    this->RepresentationProxy = 0;
  }
};

//-----------------------------------------------------------------------------
pqMIPDisplayPanelDecorator::pqMIPDisplayPanelDecorator(
  pqDisplayPanel* _panel):Superclass(_panel)
{
  pqDisplayProxyEditor *panel = qobject_cast<pqDisplayProxyEditor*> (_panel);
  pqRepresentation     *repr  = panel->getRepresentation();
  vtkSMProxy      *reprProxy  = (repr) ? repr->getProxy() : NULL;
  this->Internals             = NULL;

  //
  // If the representation doesn't have this property, then it's not our MIP representation
  //
  vtkSMProperty* prop = reprProxy->GetProperty("ActiveScalars");
  if (!prop)  {
    return;
  }

  QWidget* wid = new QWidget(panel);
  this->Internals = new pqInternals(this);
  this->Internals->Frame = wid;
  this->Internals->setupUi(wid);
  QVBoxLayout* l = qobject_cast<QVBoxLayout*>(panel->layout());
  l->addWidget(wid);
  //
  this->Internals->RepresentationProxy = vtkSMPVRepresentationProxy::SafeDownCast(reprProxy);
  this->Internals->PipelineRepresentation = qobject_cast<pqPipelineRepresentation*>(repr);

  //
  // Type
  //
  prop = reprProxy->GetProperty("TypeScalars");
  pqFieldSelectionAdaptor* adaptor= new pqFieldSelectionAdaptor(
    this->Internals->MIPTypeArray, prop);
  this->Internals->Links.addPropertyLink(
    adaptor, "attributeMode", SIGNAL(selectionChanged()),
    reprProxy, prop, 0);
  this->Internals->Links.addPropertyLink(
    adaptor, "scalar", SIGNAL(selectionChanged()),
    reprProxy, prop, 1);
  prop->UpdateDependentDomains();

  //
  // Active
  //
  prop = reprProxy->GetProperty("ActiveScalars");
  adaptor = new pqFieldSelectionAdaptor(
    this->Internals->MIPActiveArray, prop);
  this->Internals->Links.addPropertyLink(
    adaptor, "attributeMode", SIGNAL(selectionChanged()),
    reprProxy, prop, 0);
  this->Internals->Links.addPropertyLink(
    adaptor, "scalar", SIGNAL(selectionChanged()),
    reprProxy, prop, 1);
  prop->UpdateDependentDomains();

  //
  // Colour scalars display control
  //
//  this->Internals->ColorBy->setPropertyArrayName("ColorArrayName");
//  this->Internals->ColorBy->setPropertyArrayComponent("ColorAttributeType");
//  this->Internals->ColorBy->setRepresentation(this->Internals->PipelineRepresentation);
  //
  // 
  //
  this->Internals->VTKConnect->Connect(
      this->Internals->RepresentationProxy->GetProperty("Representation"),
      vtkCommand::ModifiedEvent, this, SLOT(representationTypeChanged()));

  this->Internals->Links.addPropertyLink(
    this->Internals->MIPTypeActive, "checked", SIGNAL(toggled(bool)),
    reprProxy, reprProxy->GetProperty("TypeActive"));

  //
  //
  //
  this->setupGUIConnections();
}
//-----------------------------------------------------------------------------
pqMIPDisplayPanelDecorator::~pqMIPDisplayPanelDecorator()
{
  delete this->Internals;
  this->Internals = 0;
}
//-----------------------------------------------------------------------------
void pqMIPDisplayPanelDecorator::setupGUIConnections()
{
  QObject::connect(this->Internals->MIPEditColorMapButton, SIGNAL(clicked()), this,
      SLOT(EditColour()), Qt::QueuedConnection);

  QObject::connect(this->Internals->MIPRepaintButton, SIGNAL(clicked()), this,
      SLOT(RepaintClicked()), Qt::QueuedConnection);

  QObject::connect(this->Internals->MIPActiveParticleType, SIGNAL(valueChanged(int)), this,
      SLOT(ActiveParticleTypeChanged(int)), Qt::QueuedConnection);

  QObject::connect(this->Internals->MIPTypeArray, SIGNAL(currentIndexChanged(int)), this,
      SLOT(UpdateParticleTypes()), Qt::QueuedConnection);
  
}
//-----------------------------------------------------------------------------
void pqMIPDisplayPanelDecorator::setRepresentation(
    pqPipelineRepresentation* repr)
{
  this->Internals->PipelineRepresentation = repr;
}
//-----------------------------------------------------------------------------
void pqMIPDisplayPanelDecorator::representationTypeChanged()
{
  if (this->Internals) {
    const char* reprType = vtkSMPropertyHelper
        ( this->Internals->RepresentationProxy, "Representation" ).GetAsString();
    if ( strcmp(  reprType, "MIP (particles)"  ) == 0 ) {
      this->Internals->Frame->setEnabled(true);
      vtkSMPropertyHelper(this->Internals->RepresentationProxy,
        "InterpolateScalarsBeforeMapping").Set(0);
      this->Internals->RepresentationProxy->UpdateVTKObjects();
    }
    else {
      this->Internals->Frame->setEnabled(false);
    }
  }
}
//-----------------------------------------------------------------------------
void pqMIPDisplayPanelDecorator::UpdateParticleTypes()
{
  vtkPVDataInformation *dataInfo = 
    this->Internals->PipelineRepresentation->getInputDataInformation();
  vtkPVDataInformation* geomInfo = 
    this->Internals->RepresentationProxy->GetRepresentedDataInformation();
  vtkPVDataSetAttributesInformation *pointInfo = 
    dataInfo->GetPointDataInformation();
  vtkPVArrayInformation *arrayInfo = pointInfo->GetArrayInformation(
    this->Internals->MIPTypeArray->currentText().toAscii().data());
  if (!arrayInfo) return;
  //
  QString valstr;
  double *range = arrayInfo->GetComponentRange(0);
  int ntypes = 1+static_cast<int>(range[1]);
  if (ntypes>9) {
    ntypes = 9;
    valstr = "(Error) Clamped to 9";
  }
  else {
    valstr = QString::number(ntypes);
  }
  this->Internals->MIPActiveParticleType->setMaximum(ntypes-1);
  this->Internals->typeslabel->setText(valstr);
}
//-----------------------------------------------------------------------------
void pqMIPDisplayPanelDecorator::ActiveParticleTypeChanged(int v)
{
  this->Internals->TableIndex = v;
  //
  vtkSMProperty* SettingsProperty = this->Internals->RepresentationProxy->GetProperty("ActiveParticleSettings");
  this->Internals->RepresentationProxy->UpdatePropertyInformation(SettingsProperty);
  QList<QVariant> ActiveParticleSettings = pqSMAdaptor::getMultipleElementProperty(SettingsProperty);
  //
  int ptype = ActiveParticleSettings.at(0).toString().toInt();
  if (ptype==this->Internals->MIPActiveParticleType->value()) {
    //
    bool active = ActiveParticleSettings.at(5).toString().toInt();
    this->Internals->MIPTypeActive->setChecked(active);
  }  
  for (int i=0; i<ActiveParticleSettings.size(); i++) {
    std::cout << ActiveParticleSettings.at(i).toString().toAscii().data() << std::endl;
  }
}
//-----------------------------------------------------------------------------
void pqMIPDisplayPanelDecorator::EditColour()
{
  pqApplicationCore       *core = pqApplicationCore::instance();
  pqLookupTableManager *lut_mgr = core->getLookupTableManager();
  pqScalarsToColors      *pqlut = this->Internals->PipelineRepresentation->getLookupTable();
  vtkSMProxy               *lut = (pqlut)? pqlut->getProxy() : 0;

//  pqSplotchColorScaleEditor editor(pqCoreUtilities::mainWidget());
//  editor.setActiveColorTable(pqlut);
//  editor.setRepresentation(this->Internals->PipelineRepresentation);
//  editor.exec();
}
//-----------------------------------------------------------------------------
void pqMIPDisplayPanelDecorator::RepaintClicked()
{
  if (this->Internals->PipelineRepresentation) {
    this->Internals->PipelineRepresentation->renderViewEventually();
  }
}
