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
#include "pqComboBoxDomain.h"

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
    this->PipelineRepresentation = 0;
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

  QWidget* wid = new QWidget(panel);
  this->Internals = new pqInternals(this);
  this->Internals->Frame = wid;
  this->Internals->setupUi(wid);
  QVBoxLayout* l = qobject_cast<QVBoxLayout*>(panel->layout());
  l->addWidget(wid);
  //
  this->setRepresentation(
    static_cast<pqPipelineRepresentation*> (panel->getRepresentation()));
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
void pqMIPDisplayPanelDecorator::setRepresentation(pqPipelineRepresentation* repr)
{
  if (this->Internals->PipelineRepresentation == repr) {
    return;
  }

  if (this->Internals->PipelineRepresentation) {
    // break all old links.
    this->Internals->Links.removeAllPropertyLinks();
  }

  this->Internals->PipelineRepresentation = repr;
  if (!repr) {
//    this->Internals->TransferFunctionDialog->hide();
    return;
  }

  vtkSMProperty* prop;
  vtkSMProxy *reprProxy  = repr->getProxy();
  this->Internals->RepresentationProxy = vtkSMPVRepresentationProxy::SafeDownCast(reprProxy);
  reprProxy->GetProperty("Input")->UpdateDependentDomains();

  //
  // Field array controls
  //
  // Type
  //
  prop = reprProxy->GetProperty("MIPTypeScalars");
  // adaptor from combo to property
  pqSignalAdaptorComboBox *adaptor = new pqSignalAdaptorComboBox(this->Internals->MIPTypeArray);
  // domain to control the combo contents
  new pqComboBoxDomain(this->Internals->MIPTypeArray, prop, "array_list");
  // link gui changes to property and vice versa
  this->Internals->Links.addPropertyLink(adaptor, "currentText", SIGNAL(currentTextChanged(const QString&)), reprProxy, prop);
  prop->UpdateDependentDomains();

  //
  // Active
  //
  prop = reprProxy->GetProperty("MIPActiveScalars");
  // adaptor from combo to property
  adaptor = new pqSignalAdaptorComboBox(this->Internals->MIPActiveArray);
  // domain to control the combo contents
  new pqComboBoxDomain(this->Internals->MIPActiveArray, prop, "array_list");
  // link gui changes to property and vice versa
  this->Internals->Links.addPropertyLink(adaptor, "currentText", SIGNAL(currentTextChanged(const QString&)), reprProxy, prop);
  prop->UpdateDependentDomains();

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
  vtkSMProperty* SettingsProperty = this->Internals->RepresentationProxy->GetProperty("MIPActiveParticleSettings");
  this->Internals->RepresentationProxy->UpdatePropertyInformation(SettingsProperty);
  QList<QVariant> ActiveParticleSettings = pqSMAdaptor::getMultipleElementProperty(SettingsProperty);
  //
  int ptype = ActiveParticleSettings.at(0).toString().toInt();
  if (ptype==this->Internals->MIPActiveParticleType->value()) {
    //
    bool active = ActiveParticleSettings.at(1).toString().toInt();
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
